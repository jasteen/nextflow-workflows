//sort_and_fastq.nf

ch_inputFiles = Channel.fromPath("./*.bam").map{file -> tuple(file.name.take(file.name.lastIndexOf('.')), file)}

process namesort {

    input:
        set baseName, file(bam) from ch_inputFiles
    output:
        set baseName, file("${baseName}.qsort.bam") into ch_sortedBams
            
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        4
    module      'samtools'
    memory      globalMemoryM
    time        '3h'
    queue       globalQueueL

    script:
    """
    samtools sort -n -@${task.cpus} ${bam} ${baseName}.qsort.bam
    """
}

process makeFastq {

    input:
        set baseName, file(bam) from ch_Bams
    output:
        set baseName, file("${baseName}.R1.fastq"), file("${baseName}.R2.fastq") into ch_sortedBams
            
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'bedtools'
    memory      globalMemoryM
    time        '3h'
    queue       globalQueueL

    script:
    """
    bedtools bamtofastq -i ${bam} \
                      -fq ${baseName}.R1.fastq \
                      -fq2 ${baseName}.R2.fastq
    """
}


ch_CollectedFastqs = ch_sortedBams.collect()

process gzipFastq {

    publishDir path: './processed', mode: 'copy'

    input:
        file(fastq) from ch_CollectedFastqs
    output:
        file("${fastq}.gz") into ch_zippedFastqs
            
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        '3h'
    queue       globalQueueL

    script:
    """
    gzip fastq
    """
}
