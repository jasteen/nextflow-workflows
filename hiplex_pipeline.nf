#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/home/jste0021/vh83_scratch/projects/COMPLEXO_GENESIS_NF/fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files

vardictBed       = file("/projects/vh83/reference/hiplex_complexo/COMPLEXO22_621716iii.8col.bed")
intervalFile     = file("/projects/vh83/reference/hiplex_complexo/COMPLEXO22_621716iii_final_b37.bed")

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")

AF_THR           = 0.1

//Annotation resources
dbsnp_b37       = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")
other_vep       = file("/usr/local/vep/90/ensembl-vep/cache")
vep_brcaex      = file("/projects/vh83/reference/annotation_databases/BRCA-Exchange/BRCA-exchange_accessed-180118/BRCA-exchange_accessed-180118.sort.vcf.gz")
vep_gnomad      = file("/projects/vh83/reference/annotation_databases/gnomAD/gnomad.exomes.r2.0.2.sites/gnomad.exomes.r2.0.2.sites.vcf.gz")
vep_revel       = file("/projects/vh83/reference/annotation_databases/REVEL/REVEL-030616/revel_all_chromosomes.vcf.gz")
vep_maxentscan  = file("/projects/vh83/reference/annotation_databases/MaxEntScan/MaxEntScan_accessed-240118")
vep_exac        = file("/projects/vh83/reference/annotation_databases/ExAC/ExAC_nonTCGA.r0.3.1/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz")
vep_dbnsfp      = file("/projects/vh83/reference/annotation_databases/dbNSFP/dbNSFPv2.9.3-VEP/dbNSFP-2.9.3.gz")
vep_dbscsnv     = file("/projects/vh83/reference/annotation_databases/dbscSNV/dbscSNV1.0-VEP/dbscSNV.txt.gz")
vep_cadd        = file("/projects/vh83/reference/annotation_databases/CADD/CADD-v1.3/1000G_phase3.tsv.gz")

// Tools
picardJar      = '~/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'


// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	  = 12
vepCores          = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Creating channel from input directory
ch_inputFiles = Channel.fromFilePairs("$inputDirectory/*.R{1,2}.fastq.gz")


process align_bwa {

    input:
        set baseName, file(fastqs) from ch_inputFiles
    output:
        set baseName, file("${baseName}.hq.sorted.bam"), file("${baseName}.hq.sorted.bam.bai") into ch_mappedBams

    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module      samtoolsModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL

    """
    bwa mem -M -t $bwaCores -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" $ref ${fastqs[0]} ${fastqs[1]}  \
        | samtools view -u -h -q 1 -f 2 -F 4 -F 8 -F 256 - \
        | samtools sort -@ $bwaCores -o "${baseName}.hq.sorted.bam"
    samtools index "${baseName}.hq.sorted.bam" "${baseName}.hq.sorted.bam.bai"
    """
}

process run_vardict {

    input:
        set baseName, file(bam), file(bai) from ch_mappedBams                
    output: 
        set baseName, file("${baseName}.vcf") into ch_vardictVCFs           
    
    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    module      bwaModule
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL 

    """
    export PATH=/home/jste0021/scripts/git_controlled/VarDict:$PATH
    vardict -G $ref -f $AF_THR -N $baseName -b $bam -c 1 -S 2 -E 3 -g 4 $vardictBed | \
        teststrandbias.R | \
        var2vcf_valid_b37_chrnames.pl -N $baseName -E -f $AF_THR > "${baseName}.vcf"
    """
}


//collectedVCF = Channel.from(vardictVCFS).collect().println{ "./bam_out/file.list"}


process collect_vcf_result {

      publishDir './variants_merged_out/', mode: 'copy'

      input:
      file all_vcf from ch_vardictVCFs.collect()
      file refFai

      output:
      file "merged.vcf" into ch_mergedVCF

      when:
      !all_vcf.empty

      shell:
      '''
      mkdir VCF
      mv *.vcf VCF
      # Extract the header from the first VCF
      sed '/^#CHROM/q' VCF/!{all_vcf[0]} > header.txt
      # Add contigs in the VCF header
      cat !{refFai} | cut -f1,2 | sed -e 's/^/##contig=<ID=/' -e 's/[	 ][	 ]*/,length=/' -e 's/$/>/' > contigs.txt
      sed -i '/##reference=.*/ r contigs.txt' header.txt
      # Add version numbers in the VCF header
      sed -i '/##source=.*/ r versions.txt' header.txt
      # Check if sort command allows sorting in natural order (chr1 chr2 chr10 instead of chr1 chr10 chr2)
      if [ `sort --help | grep -c 'version-sort' ` == 0 ]
      then
          sort_ops="-k1,1d"
      else
          sort_ops="-k1,1V"
      fi
      # Add all VCF contents and sort
      grep --no-filename -v '^#' VCF/*.vcf | LC_ALL=C sort -t '	' $sort_ops -k2,2n >> header.txt
      mv header.txt merged.vcf
      '''
}










//process mergeVCFS {
//    input:
//        file(vcfs) from vardictVCFs.collect().toList()
//    output:
//        file("merged.vcf.gz") into mergedVCF
//
//    publishDir path: './variants_merged_out', mode: 'copy'
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    module      bwaModule
//    memory      globalMemoryM
//    time        globalTimeM
//    queue       globalQueueL
//
//    shell:
//    '''
//    # Extract the header from the first VCF
//    grep '^#' !{all_vcf[0]} > !{out_annotated_vcf}
//    # Add version numbers in the VCF header just after fileformat line
//    grep -h -v '^#' *.vcf >> !{out_annotated_vcf}
//
//
//}
//


//cript:
//    """
//    echo $vardictVCFs.join("\n") >> subvcf.list
//    java -Xmx16g -jar ~/picard.jar mergeVcfs \
//          -I subvcf.list \
//          -O merged.vcf.gz
//    """
//}


process vt_decompose_normalise {
        
    input:
        set baseName, file(vcf) from ch_mergedVCF
    output:
        set baseName, file("merged.vt.vcf.gz") into ch_vtDecomposeVCF

    publishDir path: './variants_merged_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'vt'
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL



    """
    vt decompose -s $vcf | vt normalize -r $ref -o merged.vt.vcf.gz -
    """
}

process apply_vep {

    input:
        set baseName, file(vcf) from ch_vtDecomposeVCF
    output:
        set baseName, file("merged.vt.vep.vcf.gz") into ch_vepVCF

    publishDir path: './variants_merged', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        vepCores
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --custom $vep_brca,brcaex,vcf,exact,0,Clinical_significance_ENIGMA,\
                      Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,\
                      Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
                      --custom $vep_gnomad,gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
                      --custom $vep_revel,RVL,vcf,exact,0,REVEL_SCORE \
                      --plugin MaxEntScan,$vep_maxentscan \
                      --plugin ExAC,$vep_exac,AC,AN \
                      --plugin dbNSFP,$vep_dbnsfp,REVEL_score,REVEL_rankscore \
                      --plugin dbscSNV,$vep_dbscsnv \
                      --plugin CADD,$vep_cadd \
                      --fork $vepCores \
                      -i $vcf \
                      -o merged.vt.vep.vcf.gz}
    """
}
