#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/scratch/vh83/projects/medha_exomes/testing_fgbio/fastqs/')
outputDir      = "/scratch/vh83/projects/medha_exomes/out"
panel_bed      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions.bed')
padded_bed     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded.bed')

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
bwaCores	  = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Creating channel from input directory
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz")
inputIndexes = Channel.fromPath("$inputDirectory/*_I2.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')), file)}

process surecallTrimmer {
    input:
        set baseName, file(fastqs) from inputFiles
    output:
        set baseName, file("${baseName}_R{1,2}.fastq.gz*") into processedInputFiles,forFastqc

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    java -Xmx4g -jar ${surecalltrimmerJar} \
       -fq1 ${fastqs[0]} -fq2 ${fastqs[1]} -xt -out_loc .  
    """
       
}


process fastQC {
    publishDir path: './fastqc_results', mode: 'copy'
    input:
        set baseName, file(fastq) from forFastqc
    output:
        set baseName, file("${baseName}*.zip"), file("${baseName}*.html") into fastQCOutput

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'fastqc'
    memory      globalMemoryS
    time        globalTimeM
    queue       globalQueueL

    """
    fastqc $fastq
    """
}

process alignBwa {
    input:
        set baseName, file(fastqs) from processedInputFiles
    output:
        set baseName, file("${baseName}.bam") into bamFiles

    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module	'samtools'
    cpus        bwaCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    set -o pipefail
    bwa mem \
        -K 100000000 -v 3 -Y -t $bwaCores \
        -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" \
        $ref ${fastqs[0]} ${fastqs[1]} | \
        samtools view -bh -o ${baseName}.bam
    """
}

process sortBWA {
    input:
        set baseName, file(bam) from bamFiles
    output:
        set baseName, file("${baseName}.sorted.bam") into sortedBams
    
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module      'samtools'
    cpus        bwaCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    samtools sort -@ $bwaCores -o ${baseName}.sorted.bam $bam
    """
}

umiMap = sortedBams.join(inputIndexes)

//locatitInput.println()

process AnnotateBamWithUmis {
    input:
        set baseName, file(bam), file(index) from umiMap
    output:
        set baseName, file("${baseName}.umi.bam") into umiMarkedBams

    publishDir path: './bam_out', mode: 'copy'    

    executor    globalExecutor
    stageInMode globalStageInMode
    module	'java'
    module      'fgbio'
    cpus        globalCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    java -Xmx62g -jar '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar' AnnotateBamWithUmis \
         -i ${bam} -f ${index} -o "${baseName}.umi.bam" 
    """

}

process markDuplicates {
    input:
        set baseName, file(bam) from umiMarkedBams 
    output:
        set baseName, file("${baseName}.umi.dedup.bam"), file("${baseName}.dedup.metrics.txt") into dedupBams 
    
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module 	'picard'
    cpus        globalCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    java -Xmx62g -jar /usr/local/picard/2.9.2/bin/picard.jar MarkDuplicates I=$bam O=${baseName}.umi.dedup.bam M="${baseName}.dedup.metrics.txt" BARCODE_TAG=RX
    """

}

process groupreadsByUmi {
    input:
        set baseName, file(bam) from dedupBams
    output:
        set baseName, file("${baseName}.umi.histogram.tsv"), file("${baseName}.umi.grouped.bam") into umiGroupedBams
    
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    module      'fgbio'
    cpus        globalCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    java -Xmx62g -jar '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar' GroupReadsByUmi \
         -i ${bam} -f "${baseName}.umi.histogram.tsv" -o "${baseName}.umi.grouped.bam" -s Adjacency -e 1 
    """

}


process indexUmiBam {
    input:
        set baseName, file(hist), file(bam) from umiGroupedBams
    output:
        set baseName, file(hist), file(bam), file("${baseName}.umi.grouped.bam.bai") into indexedBams

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'samtools'
    cpus        globalCores
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    samtools index $bam ${baseName}.umi.grouped.bam.bai
    """

}


tumor  = Channel.create()
normal = Channel.create()

indexedBams.choice(tumor, normal){ a -> a[0] =~ /_FFPE$/ ? 0 : 1 }


//create bedfile segments
//bedSegments = Channel.fromPath("$padded_bed").splitText( by: 10000, file: "paddedBed.bed")
//
//

process runVardict {
    input:
        set tbaseName, file(thist), file(tbam), file(tbai) from tumor
        set nbaseName, file(nhist), file(nbam), file(nbai) from normal
    output:
        set tbaseName, nbaseName, file("${tbaseName}.${nbaseName}.somatic.vardict.raw.tsv") into rawVardictSegments
    
    publishDir path: './vardict_out', mode: 'copy'    

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryL
    time        globalTimeL
    queue       globalQueueL

    """
    /home/jste0021/scripts/git_controlled/VarDictJava/build/install/VarDict/bin/VarDict \
        -G $ref -f 0.005 -N $tbaseName -b "${tbam}|${nbam}" -c 1 -S 2 -E 3 -g 4 $padded_bed > "${tbaseName}.${nbaseName}.somatic.vardict.raw.tsv"
    """
}

//
//process runVardict {
//    input:
//        set tbaseName, file(tbam), file(tbai) from tumor
//        set nbaseName, file(nbam), file(nbai) from normal
//    output:
//        set tbaseName, nbaseName, file("${tbaseName}.${nbaseName}.somatic.vardict.raw.tsv") into rawVardict
//    
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    cpus        1
//    memory      globalMemoryL
//    time        globalTimeL
//    queue       globalQueueL
//    
//    """
//    /home/jste0021/scripts/git_controlled/VarDictJava/build/install/VarDict/bin/VarDict \
//        -G $ref -f 0.005 -N $tbaseName -b "${tbam}|${nbam}" -c 1 -S 2 -E 3 -g 4 $padded_bed > "${tbaseName}.${nbaseName}.somatic.vardict.raw.tsv" 
//    """
//}
//
//
//
//
//
//
////rawVardictSegments
//
//
////    .flatMap { key, key2, whatever, files -> files }
////    .collect()
////    .set{temp}
//
////temp.println()
//
//
////process mergeVardicts {
////    input:
////        file(segment) from temp
////    output:
////        file("somatic.vardict.raw.merged.tsv") into mergedRawVardict
////    
////
////    """
////    for every file in $segment
////    do
////        cat \$i >> somatic.vardict.raw.merged.tsv
////    done
////
////  """
////}
//
////mergedRawVardict.println()
//
//
//
//
////    """
////    cat ${beds} > $tbasebame
////                  '~/scripts/git_controlled/VarDict/teststrandbias.R | ' \
////                  '~/scripts/git_controlled/VarDict/var2vcf_valid_b37_chrnames.pl -N {sample_name} -E -f {AF_THR} > {vcf_out}'.
////    """
////}
//
//
//
//
//
//
//
//
//
//
