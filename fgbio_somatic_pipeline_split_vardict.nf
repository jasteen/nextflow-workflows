#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/scratch/vh83/projects/medha_exomes/test_split/fastqs')
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
//i'm sure there is a better way to map the basename, but this works for the moment.
ch_inputIndexes = Channel.fromPath("$inputDirectory/*_I2.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')), file)}

//join input files and index on the baseName
ch_umiMap = ch_inputFiles.join(ch_inputIndexes)


process createUnmappedUMIBam {
    
    publishDir path: './bam_out', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2), file(I2) from ch_umiMap
    output:
        set baseName, file("${baseName}.unmapped.umi.bam") into ch_unmappedUMIbams

    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'fgbio'
    module      'java'
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    """
    java -Xmx30g -Djava.io.tmpdir=$tmp_dir -XX:+AggressiveOpts -XX:+AggressiveHeap \
        -jar $fgbioJar FastqToBam --input $R1 $R2 $I2 --output "${baseName}.unmapped.umi.bam" --read-structures +T +T +M \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}


process markAdaptors {

    publishDir path: './bam_out', mode: 'copy'

    input:
        set baseName, file(bam) from ch_unmappedUMIbams
    output:
        set baseName, file("${baseName}.unmapped.umi.marked.bam"),
                      file("${baseName}.unmapped.umi.marked_metrics.tsv") into ch_markedUMIbams

    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'java'
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

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
        set baseName, file("${baseName}.piped.bam") into ch_pipedBams

    publishDir path: './bam_out', mode: 'copy'

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

    """
    set -o pipefail
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ '/dev/stdout' -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 \
        -INTERLEAVE true -NON_PF true -TMP_DIR "$tmp_dir" | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.piped.bam" -R "$ref" -CREATE_INDEX true -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}


process groupreadsByUmi {
    input:
        set baseName, file(bam) from ch_pipedBams
    output:
        set baseName, file("${baseName}.piped.grouped.histogram.tsv"), file("${baseName}.piped.grouped.bam") into ch_umiGroupedBams
    
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    """
    java -Xmx6g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar GroupReadsByUmi \
         -i ${bam} -f "${baseName}.piped.grouped.histogram.tsv" -o "${baseName}.piped.grouped.bam" -s Adjacency -e 1 
    """

}


process generateConsensusReads {
    input:
        set baseName, file(hist), file(bam) from ch_umiGroupedBams
    output:
        set baseName, file("${baseName}.consensus.unmapped.bam") into ch_unmappedConsensusBams
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL


    """
    java -Xmx6g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar CallMolecularConsensusReads \
        --input $bam --output ${baseName}.consensus.unmapped.bam \
        --error-rate-post-umi 30 --min-reads 1
    """
}

//process filterConsensusReads {
//    input:
//        set baseName, file(bam) from unmappedConsensusBams 
//    output:
//        set baseName, file("${baseName}.consensus.filtered.unmapped.bam") into filteredConsensusBams
//    publishDir path: './bam_out', mode: 'copy'
//
//    executor    globalExecutor
//    stageInMode globalStageInMode
//    module      'java'
//    module      'fgbio'
//    cpus        globalCores
//    memory      globalMemoryM
//    time        globalTimeM
//    queue       globalQueueL
//
//
//    """
//    java -Xmx6g -jar '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar' FilterConsensusReads \
//        --input $bam --output ${baseName}.consensus.filtered.unmapped.bam \
//        --ref $ref --reverse-per-base-tags true --min-reads 1 \
//        -E 0.05 -N 40 -e 0.1 -n 0.1
//    """
//
//}


process mapConsensusReads {
     input:
        set baseName, file(bam) from ch_unmappedConsensusBams
    output:
        set baseName, file("${baseName}.consensus.aligned.bam") into ch_mappedConsensusBams
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    module 	    bwaModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL


    """
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR $tmp_dir | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.consensus.aligned.bam" -R $ref -ADD_MATE_CIGAR true \
        -SO coordinate -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """

}


process indexBam {
    input:
        set baseName, file(bam) from ch_mappedConsensusBams
    output:
        set baseName, file(bam), file("${baseName}.consensus.aligned.bam.bai") into ch_indexedConsensusBams
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'samtools'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL


    """
    samtools index $bam ${baseName}.consensus.aligned.bam.bai
    """

}

//I cant think of a better workflow to make sure Tumor and Normal are processed in order
//create channel for normal and tumor
ch_tumor  = Channel.create()
ch_normal = Channel.create()

//split single bam channel into tumor and normal **CURRENTLY RELIES ON SAMPLE_[FFPE|NORMAL] naming scheme
ch_indexedConsensusBams.choice(ch_tumor, ch_normal){ a -> a[0] =~ /_FFPE$/ ? 0 : 1 }

//split SAMPLE from FFPE|NORMAL so channels can be joined by sample
ch_normalSplit = ch_normal.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}
ch_tumorSplit = ch_tumor.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}

//merge tumor and normal back together by sample number. 
ch_vardictInput = ch_tumorSplit.join(ch_normalSplit)

//TODO:  write split up vardict process by bed region to make it quicker.  will require writing tsv files 
//       and then merging before making VCF file.

//create bedfile segments
bedSegments = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")
//bedSegments = Channel.fromPath("$padded_bed").splitText( by: 50000)


process runVardict {
    input:
        set sample, ttype, file(tbam), file(tbai), ntype, file(nbam), file(nbai) from ch_vardictInput
        each file(segment) from bedSegments   
    output:
        set sample, tbam, nbam, file("${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv") into ch_rawVardictSegments
    
    publishDir path: './bam_out', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    """
    export PATH=/home/jste0021/scripts/git_controlled/VarDict:/home/jste0021/scripts/git_controlled/VarDictJava/build/install/VarDict/bin:$PATH
    VarDict -G ${ref} -f 0.01 -N "${tbam}|${nbam}" \
        -b "${tbam}|${nbam}" -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv"
    """ 

}

//ch_rawVardictSegments.println()


ch_collatedSegments = ch_rawVardictSegments.map{ sample, tbam, nbam, segment -> [sample, tbam, nbam, segment] }.groupTuple(by: [0,1,2])

ch_collatedSegments.println()

/*
process catSegments {
    input: set sample, files("*.tsv") from ch_collatedSegments.collect()
    output: set "$sample.name.toString().tokenize('_').get(0)", "$sample.name.toString().tokenize('_').get(1)", "$sample.name.toString().tokenize('_').get(2)", file("${sample}.collated.tsv") into ch_rawVardict

    publishDir path: './bam_out', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    """
    cat $files > ${sample}.collated.vardict.tsv
    """

}

process makeVCF {
    input:
        set sample, file(tsv) from ch_rawVardict
    output:
        set sample, tbam, nbam, file("${sample}.somatic.vardict.vcf") into ch_outputVCF
    
    publishDir path: './bam_out', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'R/3.5.1'
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    """  
    cat $tsv | /home/jste0021/scripts/git_controlled/VarDict/testsomatic.R | \
    var2vcf_paired.pl -N "${tbam}|${nbam}" -f 0.01 > "${sample}.:wqsomatic.vardict.vcf"
    """
}
*/
