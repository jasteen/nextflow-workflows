#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
panel_int      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions_b37.interval_list')
padded_int     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded_b37.interval_list')
panel_bed      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions_b37.bed')
padded_bed     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded_b37.bed')
tmp_dir        = file('/scratch/vh83/tmp/')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

//VEP
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
picardJar          = '~/picard.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
gatkModule         = 'gatk/4.0.11.0' 
rModule            = 'R/3.5.1'          
fgbioJar           = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	      = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Creating channel from input directory
//create channel flat because we want to join it later, and the tuple makes that more annoying than I want it to be
ch_inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz", flat: true)
//i'm sure there is a better way to map the basename, but this works for the moment.
ch_inputIndexes = Channel.fromPath("$inputDirectory/*_I2.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')), file)}

//join input files and index on the baseName
ch_umiMap = ch_inputFiles.join(ch_inputIndexes)

process createUnmappedUMIBam {
    
    publishDir path: './output/intermediate', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2), file(I2) from ch_umiMap
    output:
        set baseName, file("${baseName}.unmapped.umi.bam") into ch_unmappedUMIbams

    publishDir path: './output/intermediate', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'fgbio'
    module      'java'
    memory      globalMemoryM
    time        '6h'
    queue       globalQueueL

    script:
    """
    java -Xmx30g -Djava.io.tmpdir=$tmp_dir -XX:+AggressiveOpts -XX:+AggressiveHeap \
        -jar $fgbioJar FastqToBam --input $R1 $R2 $I2 --output "${baseName}.unmapped.umi.bam" --read-structures +T +T +M \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}

process markAdaptors {

    publishDir path: './output/intermediate', mode: 'copy'

    input:
        set baseName, file(bam) from ch_unmappedUMIbams
    output:
        set baseName, file("${baseName}.unmapped.umi.marked.bam"),
                      file("${baseName}.unmapped.umi.marked_metrics.tsv") into ch_markedUMIbams

    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'java'
    memory      globalMemoryM
    time        '3h'
    queue       globalQueueL

    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx30g -jar $picardJar MarkIlluminaAdapters \
        -INPUT $bam \
        -OUTPUT "${baseName}.unmapped.umi.marked.bam" \
        -METRICS "${baseName}.unmapped.umi.marked_metrics.tsv"
    """
}

process alignBwa {
    input:
        set baseName, file(bam), file(metrics) from ch_markedUMIbams
    output:
        set baseName, file("${baseName}.aligned.bam") into ch_pipedBams, ch_mappedNoUMI, ch_forMetrics1

    publishDir path: './output/intermediate', mode: 'copy'

    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module	    'samtools'
    module      'picard'
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:
    """
    set -o pipefail
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ '/dev/stdout' -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 \
        -INTERLEAVE true -NON_PF true -TMP_DIR "$tmp_dir" | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.aligned.bam" -R "$ref" -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process indexPreUmiBam {
    input:
        set baseName, file(bam) from ch_mappedNoUMI
    output:
        set baseName, file(bam), file("${baseName}.aligned.bam.bai") into ch_indexedMappedNoUMI
    publishDir path: './output/intermediate', mode: 'copy'

    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    module      'samtools'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:
    """
    samtools index $bam ${baseName}.aligned.bam.bai
    """

}

ch_tumorPREUMI  = Channel.create()
ch_normalPREUMI = Channel.create()

//split single bam channel into tumor and normal **CURRENTLY RELIES ON "SAMPLE_[FFPE|NORMAL]" naming scheme
ch_indexedMappedNoUMI.choice(ch_tumorPREUMI, ch_normalPREUMI){ a -> a[0] =~ /FFPE$/ ? 0 : 1 }

//split SAMPLE from FFPE|NORMAL so channels can be joined by sample
ch_normalSplitPREUMI = ch_normalPREUMI.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}
ch_tumorSplitPREUMI = ch_tumorPREUMI.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}

//merge tumor and normal back together by sample number. 
ch_tumorNormalPairsPREUMI = ch_tumorSplitPREUMI.join(ch_normalSplitPREUMI)


ch_bedSegments2 = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")
ch_vardictPREUMI= ch_tumorNormalPairsPREUMI.combine(ch_bedSegments2)

process runVardictPREUMI {
    input:
        set sample, ttype, file(tbam), file(tbai), ntype, file(nbam), file(nbai), file(segment) from ch_vardictPREUMI
    output:
        set sample, file(tbam), file(nbam), file("${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv") into ch_rawVardictSegmentsPREUMI

    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        6
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    errorStrategy 'retry'
    maxRetries  3

    script:
    """
    export PATH=/home/jste0021/scripts/VarDict-1.5.8/bin/:$PATH
    VarDict -G ${ref} -f 0.01 -N "${tbam}|${nbam}" \
        -b "${tbam}|${nbam}" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv"
    """ 

}


ch_collatedSegmentsPREUMI = ch_rawVardictSegmentsPREUMI.map{ sample, tbam, nbam, segment -> [sample, tbam.name, nbam.name, segment]}.groupTuple(by: [0,1,2])


process catSegmentsPREUMI {
    echo true
    input: 
        set sample, tbam, nbam, file(tsv) from ch_collatedSegmentsPREUMI
    output: 
        set sample, tbam, nbam, file("${sample}.collated.vardict.tsv") into ch_rawVardictPREUMI

    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    script:
    
    myfiles = tsv.collect().join(' ')

    """
    cat ${myfiles} > ${sample}.collated.vardict.tsv
    """

}

process makeVCFPREUMI {
    input:
        set sample, tbam, nbam, file(tsv) from ch_rawVardictPREUMI
    output:
        set sample, file("${sample}.somatic.vardict.vcf") into ch_outputVCFPREUMI
    
    publishDir path: './output/preUMI/intermediate', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:
    """  
    module purge
    module load R/3.5.1
    cat $tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/testsomatic.R | \
        /home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_paired.pl -N "${tbam}|${nbam}"  \
        -f 0.01 > "${sample}.somatic.vardict.vcf"
    """
}

process reheaderPREUMIVCF {
    input:
        set sample, file(vcf) from ch_outputVCFPREUMI
    output:
        set sample, file("*.vcf.gz") into ch_reheaderVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    module      'bcftools/1.8'

    script:
    """
    bcftools annotate -h ~/vh83/reference/genomes/b37/vcf_contig_header_lines.txt -O v ${vcf} | \
        bcftools sort -o ${sample}.varditc.sorted.vcf.gz -O z -
    """

}

process sortVCFSPREUMI {

    input:
        set baseName, file(vcf) from ch_reheaderVCFPREUMI
    output:
        set baseName, file("${baseName}.reheader.sorted.vcf.gz") into ch_sortedVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'                                    
    
    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    script:
    """
    bcftools sort -o "${baseName}.reheader.sorted.vcf.gz" -O z ${vcf}
    """
}

process indexVCFSPREUMI {
    input:
        set baseName, file(vcf) from ch_sortedVCFPREUMI
    output:
        set baseName, file(vcf), file("${baseName}.reheader.sorted.vcf.gz.tbi") into ch_indexedVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'                                    
    
    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    script:
    """
    bcftools index -f --tbi ${vcf} -o ${baseName}.reheader.sorted.vcf.gz.tbi
    """
}

process vt_decompose_normalisePREUMI {
        
    input:
        set baseName, file(vcf), file(tbi) from ch_indexedVCFPREUMI
    output:
        set baseName, file("${baseName}.reheader.sorted.vt.vcf.gz") into ch_vtDecomposeVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'vt'
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL

    """
    vt decompose -s $vcf | vt normalize -r $ref -o "${baseName}.reheader.sorted.vt.vcf.gz" -
    """
}

process apply_vepPREUMI {

    input:
        set baseName, file(vcf) from ch_vtDecomposeVCFPREUMI
    output:
        set baseName, file("${baseName}.reheader.sorted.vt.vep.vcf") into ch_vepVCFPREUMI

    publishDir path: './output/preUMI/somatic_annotated', mode: 'copy'

    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        12
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL
    module      'vep/90'

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --custom $vep_brcaex,brcaex,vcf,exact,0,Clinical_significance_ENIGMA,\
                      Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,\
                      Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
                      --custom $vep_gnomad,gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
                      --custom $vep_revel,RVL,vcf,exact,0,REVEL_SCORE \
                      --plugin MaxEntScan,$vep_maxentscan \
                      --plugin ExAC,$vep_exac,AC,AN \
                      --plugin dbNSFP,$vep_dbnsfp,REVEL_score,REVEL_rankscore \
                      --plugin dbscSNV,$vep_dbscsnv \
                      --plugin CADD,$vep_cadd \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o "${baseName}.reheader.sorted.vt.vep.vcf"
    """
}
