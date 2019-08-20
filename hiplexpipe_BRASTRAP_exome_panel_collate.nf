#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files

vardictBed       = file("/projects/vh83/reference/brastrap_specific/vardict/BRASTRAP_721717_8column.bed")
intervalFile     = file("/projects/vh83/reference/brastrap_specific/BRA-STRAP_621717_100.final.roverfile_g37.numsort.sorted.bed")
restrictedBed    = file("/projects/vh83/reference/brastrap_specific/BRA-STRAP_coding_regions_targeted_sort.bed")
primer_bedpe_file= file("/projects/vh83/reference/prostrap/final_prostrap_b37_bedpe_bamclipper.txt")

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")
genome_file      = file("/projects/vh83/reference/genomes/b37/accessory_files/human_g1k_v37_decoy_GenomeFile.txt")
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

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
bamclipper_exe = '/projects/vh83/local_software/bamclipper/bamclipper.sh'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	      = 4
vepCores          = 4
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '15m'
globalTimeL       = '15m'
globalQueueS      = 'short'
globalQueueL      = 'short'

// Creating channel from input directory
ch_inputFiles = Channel.fromFilePairs("${inputDirectory}/*_R{1,2}_001.fastq.gz")


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
    bwa mem -M -t ${task.cpus} -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" $ref ${fastqs[0]} ${fastqs[1]}  \
        | samtools view -u -h -q 1 -f 2 -F 4 -F 8 -F 256 - \
        | samtools sort -@ $bwaCores -o "${baseName}.hq.sorted.bam"
    samtools index "${baseName}.hq.sorted.bam" "${baseName}.hq.sorted.bam.bai"
    """
}

ch_mappedBams.into{ch_mappedBam1;ch_mappedBam2;ch_mappedBam3;ch_mappedBam4;ch_mappedBam5;ch_forBamClipper}

process run_bamClipper {
    input:
        set baseName, file(bam), file(bai) from ch_forBamClipper               
    output: 
        set baseName, file("${baseName}.hq.sorted.primerclipped.bam"), file("${baseName}.hq.sorted.primerclipped.bam.bai") into ch_forperBase          
    
    publishDir path: './bamclipper', mode: 'copy'                                    
    
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL 
    module      'gnuparallel/20190122'
    module      'samtools'

    """
    ${bamclipper_exe} -b ${bam} -p ${primer_bedpe_file} -n ${task.cpus}
    """
}

process generatePerbaseMetrics {
    echo true
    input:
        set baseName, file(bam), file(index) from ch_forperBase
                 
    output: 
        set baseName, file("${baseName}.mpileup.sorted.vcf.gz"), file("${baseName}.mpileup.sorted.vcf.gz.tbi")  into ch_mpileupOUT                                    
    
    
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryL 
    time        globalTimeM
    queue       globalQueueL 
    module      'samtools'
    cpus        '2'
    module      'bcftools'
    

    
    """
    bcftools mpileup --threads ${task.cpus} -Ou -d 250 -B -R ${restrictedBed} \
        -a "FORMAT/DP" -f ${ref} ${bam} | bcftools sort -o "${baseName}.mpileup.sorted.vcf.gz" -O z ${vcf}
    bcftools index -f --tbi "${baseName}.sorted.vcf.gz"
    """
}


ch_mpileupOUT
    .map { mytuple -> [ mytuple.collect{ it[1] }, mytuple.collect{ it[2] } ] }
    .set{ch_fucks_given}


process mergepileupVCFS {
    echo true
    publishDir './variants_merged/', mode: 'copy'
    input:
    set file(vcf), file(index) from ch_fucks_given
        
    output:
    file "full_merged_mpileup.vcf.gz" into ch_mergedVCF

    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    
    script: 
    
    """
    ls *.vcf.gz | split -l10 split_vcf_
    for i in split_vcf_*; do bcftools merge -O z -o "\${i}.mpileup_merged.vcf.gz"  \
        -l \$i; mpileup index -f --tbi \${i}.mpileup_merged.vcf.gz; done
    ls "*.mpileup_merged.vcf.gz" > temp.txt
    bcftools merge -l temp.txt -O z -o "full_merged_mpileup.vcf.gz"
    """
}