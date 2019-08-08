//sort_and_fastq.nf

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'
picardJar          = '~/picard.jar'


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
    samtools sort -n -@${task.cpus} ${bam} -o ${baseName}.qsort.bam
    """
}

process makeFastq {

    publishDir path: './processed', mode: 'copy' 
    input:
        set baseName, file(bam) from ch_sortedBams
    output:
        file("*.fastq.gz") into ch_fastqs
            
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
    java -Dpicard.useLegacyParser=false -Xmx30g -jar $picardJar SamToFastq \
     -I $bam \
     -F ${baseName}_R1.fastq.gz \
     -F2 ${baseName}_R2.fastq.gz \
     """
}

