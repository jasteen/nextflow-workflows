#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/scratch/vh83/projects/medha_exomes/fastqs/')
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
globalMemoryL     = '16 GB'
globalTimeS       = '8m'
globalTimeL       = '8h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Creating channel from input directory
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz", flat:true).splitFastq(by:20000000, file:true, pe:true, compress:true)
inputIndexes = Channel.fromPath("$inputDirectory/*_I2.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')-3), file)}


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

processedInputFiles.println()

//process fastQC {
//    publishDir path: '/fastqc_results', mode: 'copy'
//    input:
//        set baseName, file(fastq) from forFastqc
//    output:
//        set baseName, file("${baseName}*.zip"), file("${baseName}*.html") into fastQCOutput
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    cpus        1
//    module      'fastqc'
//    memory      globalMemoryS
//    time        globalTimeS
//    queue       globalQueueS
//
//    """
//    fastqc $fastq
//    """
//}
//
//process alignBwa {
//    input:
//        set baseName, file(fastqs) from processedInputFiles
//    output:
//        set baseName, file("${baseName}.bam") into bamFiles
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    module      bwaModule
//    module	'samtools'
//    cpus        globalCores
//    memory      globalMemoryM
//    time        globalTimeL
//    queue       globalQueueL
//
//    """
//    set -o pipefail
//    bwa mem \
//        -K 100000000 -v 3 -Y -t $bwaCores \
//        -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" \
//        $ref ${fastqs[0]} ${fastqs[1]} | \
//        samtools view -u -h -o ${baseName}.bam | \
//        samtools sort -@ $bwaCores -o ${baseName}.bam  
//    """
//}
//
//locatitInput = bamFiles.join(inputIndexes)
//
//process runLocatit {
//    input:
//        set baseName, file(bam), file(index) from locatitInput
//    output:
//        set baseName, file("${baseName}.locatit.bam") into locatitBams
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    module	'java'
//    cpus        globalCores
//    memory      globalMemoryM
//    time        globalTimeS
//    queue       globalQueueS
//
//    """
//    java -Xmx4g -jar '/projects/vh83/local_software/agent/LocatIt_v4.0.1.jar' \
//         -i -C -U -q 25 -m 3 -c 2500 -d 0 -IB -OB -b ${panel_bed}  \
//         -o "${baseName}.locatit.bam" ${bam} ${index}
//    """
//
//}
//
//
//process sortBam {
//    input:
//        set baseName, file(bam) from locatitBams
//    output:
//        set baseName, 
//            file("${baseName}.locatit.sorted.bam"), 
//            file("${baseName}.locatit.sorted.bai") into sortedBamFiles
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
//        INPUT=$bam \
//        OUTPUT=${baseName}.locatit.sorted.bam \
//        SORT_ORDER=coordinate \
//        CREATE_INDEX=true \
//        CREATE_MD5_FILE=true \
//        MAX_RECORDS_IN_RAM=300000
//    """
//}
//
//tumor  = Channel.create()
//normal = Channel.create()
//sortedBamFiles.choice(tumor, normal){ a -> a[0] =~ /_FFPE$/ ? 0 : 1 }
//
//
//
//
////create bedfile segments
//bedSegments = Channel.fromPath("$padded_bed").splitText( by: 10000, file: "paddedBed.bed")
//
//
//process runVardict {
//    input:
//        set tbaseName, file(tbam), file(tbai) from tumor
//        set nbaseName, file(nbam), file(nbai) from normal
//        each file(segment) from bedSegments
//    output:
//        set tbaseName, nbaseName, segment, file("${tbaseName}.${nbaseName}.${segment}.somatic.vardict.raw.tsv") into rawVardictSegments
//    
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    cpus        1
//    memory      globalMemoryS
//    time        globalTimeS
//    queue       globalQueueS
//
//    """
//    /home/jste0021/scripts/git_controlled/VarDictJava/build/install/VarDict/bin/VarDict \
//        -G $ref -f 0.005 -N $tbaseName -b "${tbam}|${nbam}" -c 1 -S 2 -E 3 -g 4 $segment > "${tbaseName}.${nbaseName}.${segment}.somatic.vardict.raw.tsv"
//    """
//}
//
//
//
//rawVardictSegments
//    .flatMap { key, key2, whatever, files -> files }
//    .collect()
//    .set{temp}
//
////temp.println()
//
//
//process mergeVardicts {
//    input:
//        file(segment) from temp
//    output:
//        file("somatic.vardict.raw.merged.tsv") into mergedRawVardict
//    
//
//    """
//    for every file in $segment
//    do
//        cat \$i >> somatic.vardict.raw.merged.tsv
//    done
//
//  """
//}
//
//mergedRawVardict.println()
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
