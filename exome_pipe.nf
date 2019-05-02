#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
panel_bed      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions.bed')
padded_bed     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded.bed')
tmp_dir        = file('/scratch/vh83/tmp/')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")

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


process createUnmappedUMIBam {
    
    publishDir path: './output/intermediate', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2), from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped.bam") into ch_unmappedBams

    publishDir path: './output/intermediate', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'fgbio'
    module      'java'
    memory      globalMemoryM
    time        '3h'
    queue       globalQueueL

    script:
    """
    java -Xmx30g -Djava.io.tmpdir=$tmp_dir -XX:+AggressiveOpts -XX:+AggressiveHeap \
        -jar $fgbioJar FastqToBam --input $R1 $R2 --output "${baseName}.unmapped.bam" --read-structures +T +T \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}

process markAdaptors {

    publishDir path: './output/intermediate', mode: 'copy'

    input:
        set baseName, file(bam) from ch_unmappedBams
    output:
        set baseName, file("${baseName}.unmapped.marked.bam"),
                      file("${baseName}.unmapped.marked_metrics.tsv") into ch_markedBams


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
        -OUTPUT "${baseName}.unmapped.marked.bam" \
        -METRICS "${baseName}.unmapped.marked_metrics.tsv"
    """
}


process alignBwa {
    input:
        set baseName, file(bam), file(metrics) from ch_markedBams
    output:
        set baseName, file("${baseName}.mapped.bam") into ch_mappedBams, ch_forMetrics

    publishDir path: './output/intermediate', mode: 'copy'

    cache       'deep'
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
        -OUTPUT "${baseName}.piped.bam" -R "$ref" -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process indexBam {
    input:
        set baseName, file(bam) from ch_mappedBams
    output:
        set baseName, file(bam), file("${baseName}.mapped.bam.bai") into ch_indexedMapped
    publishDir path: './output/intermediate', mode: 'copy'

    cache       'deep'
    executor    globalExecutor
    stageInMode globalStageInMode
    module      'samtools'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:
    """
    samtools index $bam ${baseName}.mapped.bam.bai
    """

}


ch_bedSegments = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")

ch_vardict= ch_indexedMapped.combine(ch_bedSegments)

process vardictPreUMI {
    
    input:
        set baseName, file(bam), file(bai), file(segment) from ch_vardict
    output:
        set baseName, file("${baseName}.${segment}.vardict.tsv") into ch_vardictSegments

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        6
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:
    """
    export PATH=/home/jste0021/scripts/VarDict-1.5.8/bin/:$PATH
    VarDict -G ${ref} -f 0.01 -N "$baseName" \
        -b "$bam" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${baseName}.${segment}.vardict.tsv"
    """

}

ch_collatedSegment = ch_vardictSegments.map{ sample, segment -> [sample, segment]}.groupTuple(by: [0])

process catSegments {
    echo true

    input: 
        set baseName, file(tsv) from ch_collatedSegments
    output: 
        set baseName, file("${baseName}.collated.vardict.tsv") into ch_vardictCollated
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    script:
    
    myfiles = tsv.collect().join(' ')

    """
    cat ${myfiles} > ${baseName}.collated.vardict.tsv
    """

}


process makeVCF {
    input:
        set sample, file(tsv) from ch_vardictCollated
    output:
        set sample, file("${sample}.vardict.vcf") into ch_outputVCF
    
    publishDir path: './output/variants', mode: 'copy'
    
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
    cat $tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/teststrandbias.R | \
    /home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_valid.pl -N "$sample" -f 0.01 > "${sample}.vardict.vcf"
    """
}

process collectMetrics {

    input:
        set sample, file(bam) from ch_forMetrics
    output:
        set sample, file(*multiple_metrics*) into ch_metrics
    
    publishDir path: './output/metrics', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:

    """
    java -Dpicard.useLegacyParser=false -Xmx6G -jar ${picardJar} CollectMultipleMetrics \
        -I $bam \
        -O ${sample}.multiple_metrics \
        -R $ref 
    """


}
