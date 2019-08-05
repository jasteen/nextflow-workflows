//sort_and_fastq.nf

process namesort {

    publishDir path: './output/intermediate', mode: 'copy'

    input:
        set baseName, file(bam) from ch_unmappedBams
    output:
        set baseName, file("${baseName}.unmapped.marked.bam"),
                      file("${baseName}.unmapped.marked_metrics.tsv") into ch_markedBams

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
        -OUTPUT "${baseName}.unmapped.marked.bam" \
        -METRICS "${baseName}.unmapped.marked_metrics.tsv"
    """
}

