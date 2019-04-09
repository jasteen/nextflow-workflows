#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('/scratch/vh83/projects/medha_exomes/full_fg_workflow/fastqs/')
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
picardJar      = '~/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'
gatkModule     = 'gatk/4.0.11.0' 
surecalltrimmerJar = '/projects/vh83/local_software/agent/SurecallTrimmer_v4.0.1.jar'
locatitJar     = '/projects/vh83/local_software/agent/LocatIt_v4.0.1.jar'
fgbioJar       = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	  = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Creating channel from input directory
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz", flat: true)
inputIndexes = Channel.fromPath("$inputDirectory/*_I2.fastq.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('_')), file)}

umiMap = inputFiles.join(inputIndexes)


process createUnmappedUMIBam {
    publishDir path: './bam_out', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2), file(I2) from umiMap
    output:
        set baseName, file("${baseName}.unmapped.umi.bam") into unmappedUMIbams

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
    java -Xmx30g -Djava.io.tmpdir=/scratch/vh83/tmp -XX:+AggressiveOpts -XX:+AggressiveHeap -jar $fgbioJar FastqToBam --input $R1 $R2 $I2 --output "${baseName}.unmapped.umi.bam" --read-structures +T +T +M \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}


process markAdaptors {

    publishDir path: './bam_out', mode: 'copy'

    input:
        set baseName, file(bam) from unmappedUMIbams
    output:
        set baseName, file("${baseName}.unmapped.umi.marked.bam"), file("${baseName}.unmapped.umi.marked_metrics.tsv") into markedUMIbams

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
        set baseName, file(bam), file(metrics) from markedUMIbams
    output:
        set baseName, file("${baseName}.piped.bam") into pipedBams

    publishDir path: './bam_out', mode: 'copy'

    cache       'deep'
    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module	'samtools'
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
    bwa mem -M -t $bwaCores -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.piped.bam" -R "$ref" -CREATE_INDEX true -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}


process groupreadsByUmi {
    input:
        set baseName, file(bam) from pipedBams
    output:
        set baseName, file("${baseName}.piped.grouped.histogram.tsv"), file("${baseName}.piped.grouped.bam") into umiGroupedBams
    
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    """
    java -Xmx6g -Djava.io.tmpdir=/scratch/vh83/tmp -jar '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar' GroupReadsByUmi \
         -i ${bam} -f "${baseName}.piped.grouped.histogram.tsv" -o "${baseName}.piped.grouped.bam" -s Adjacency -e 1 
    """

}


process generateConsensusReads {
    input:
        set baseName, file(hist), file(bam) from umiGroupedBams
    output:
        set baseName, file("${baseName}.consensus.unmapped.bam") into unmappedConsensusBams
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL


    """
    java -Xmx6g -Djava.io.tmpdir=/scratch/vh83/tmp -jar '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar' CallMolecularConsensusReads \
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
        set baseName, file(bam) from unmappedConsensusBams
    output:
        set baseName, file("${baseName}.consensus.aligned.bam") into mappedConsensusBams
    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      'java'
    module 	bwaModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL


    """
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR $tmp_dir | \
    bwa mem -M -t $bwaCores -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.consensus.aligned.bam" -R $ref -ADD_MATE_CIGAR true \
        -SO coordinate -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """

}


process indexBam {
    input:
        set baseName, file(bam) from mappedConsensusBams
    output:
        set baseName, file(bam), file("${baseName}.consensus.aligned.bam.bai") into indexedConsensusBams
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

//I cant think of a better workflow to make sure Tumor and normal are processed in order
//create channel for normal and tumor
tumor  = Channel.create()
normal = Channel.create()

//split single bam channel into tumor and normal
indexedConsensusBams.choice(tumor, normal){ a -> a[0] =~ /_FFPE$/ ? 0 : 1 }

//fix the names to allow remerging by sample number
normalSplit = normal.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}
tumorSplit = tumor.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}

//merge tumor and normal back together by sample number. 
vardictInput = tumorSplit.join(normalSplit)


process runVardict {
    input:
        set sample, ttype, file(tbam), file(tbai), ntype, file(nbam), file(nbai) from vardictInput
    output:
        set sample, file("${sample}.${ttype}_v_${ntype}.somatic.vardict.vcf") into rawVardict
    
    publishDir path: './bam_out', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    module      'R/3.5.1'
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL
    
    """
    export PATH=/home/jste0021/scripts/git_controlled/VarDict:$PATH
    /home/jste0021/scripts/git_controlled/VarDictJava/build/install/VarDict/bin/VarDict \
        -G $ref -f 0.01 -N "${tbam}|${nbam}" -b "${tbam}|${nbam}" -c 1 -S 2 -E 3 -g 4 $padded_bed | \
    /home/jste0021/scripts/git_controlled/VarDict/testsomatic.R | \
    var2vcf_paired.pl -N "${tbam}|${nbam}" -f 0.01 > "${sample}.${ttype}_v_${ntype}.somatic.vardict.vcf"
    """
}

