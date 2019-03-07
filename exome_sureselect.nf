#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/scratch/vh83/sandbox/jason/fastqs/')
outputDir      = "/scratch/vh83/sandbox/jason/out"
panel_bed      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions.bed')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")


// Tools
picardJar      = '/usr/local/picard/2.9.2/bin/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'
gatkModule     = 'gatk/4.0.11.0' 
surecalltrimmerJar = '/projects/vh83/local_software/agent/SurecallTrimmer_v4.0.1.jar'
locatitJar     = '/projects/vh83/local_software/agent/LocatIt_v4.0.1.jar'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '16 GB'
globalTimeS       = '8m'
globalQueueS      = 'short'


// Creating channel from input directory
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,3}_*.fastq.gz")
inputIndexes = Channel.fromPath("$inputDirectory/*_R2_*.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')-3), file)}

process surecallTrimmer {
    input:
        set baseName, file(fastqs) from inputFiles
    output:
        set baseName, file("${baseName}_R{1,3}_*.fastq.gz*") into processedInputFiles

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS

    """
    java -Xmx4g -jar ${surecalltrimmerJar} \
       -fq1 ${fastqs[0]} -fq2 ${fastqs[1]} -xt -out_loc .  
    """
       
}

process alignBwa {
    input:
        set baseName, file(fastqs) from processedInputFiles
    output:
        set baseName, file("${baseName}.bam") into bamFiles

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module	'samtools'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS

    """
    set -o pipefail
    bwa mem \
        -K 100000000 -v 3 -Y -t $globalCores \
        -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" \
        $ref ${fastqs[0]} ${fastqs[1]} | \
        samtools view -b -h -o ${baseName}.bam -
    """
}

locatitInput = bamFiles.join(inputIndexes)

process runLocatit {
    input:
        set baseName, file(bam), file(index) from locatitInput
    output:
        set baseName, file("${baseName}.locatit.bam") into locatitBams

    executor    globalExecutor
    stageInMode globalStageInMode
    module	'java'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS

    """
    java -Xmx4g -jar $locatitJar \
       -U -q 25 -m 1 -d 0 -IB -OB -b ${panel_bed} \
       -o "${baseName.locatit.bam}" ${bam} ${index}   
    """

}

locatitBams.println()

//process sortBam {
//    input:
//        set baseName, markedBam from markedBamFiles
//    output:
//        set baseName,
//            file("${baseName}.marked.sorted.bam"), 
//            file("${baseName}.marked.sorted.bai") into sortedBamFiles
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    cpus        1
//    memory      globalMemoryS
//    time        globalTimeS
//    queue       globalQueueS
//
//    """
//    java -Xmx4000m -jar $picardJar SortSam \
//        INPUT=$markedBam \
//        OUTPUT=${baseName}.marked.sorted.bam \
//        SORT_ORDER=coordinate \
//        CREATE_INDEX=true \
//        CREATE_MD5_FILE=true \
//        MAX_RECORDS_IN_RAM=300000
//    """
//}
//
//
