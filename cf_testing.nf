#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./test')
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
vep_cache      = file("/projects/vh83/reference/VEP_CACHE")


// Tools
picardJar          = '~/picard.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
gatkModule         = 'gatk/4.0.11.0' 
rModule            = 'R/3.5.1'          
fgbioJar           = '/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar'
condaModule        = 'miniconda3/4.1.11-python3.5' 

// Creating channel from input directory
Channel.fromFilePairs("$inputDirectory/*_{R1,R2}.fastq.gz", size: 2, flat: true).into{ch_inputFiles;ch_forFastqc}

process runFASTQC {
    label 'genomics_3'

    input:
        set baseName, file(R1), file(R2) from ch_forFastqc
    output:
        file("*.{html,zip}") into ch_fastqcReports
    
    publishDir path: './output/metrics/fastqc', mode: 'copy'

    module      'fastqc'
        
    script:
    
    """
    fastqc -t ${task.cpus} -q $R1 $R2
    """

}

process surecallTrimmer {
    
    label 'genomics_3'

    input:
        set baseName, file(R1), file(R2) from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped_R1.fastq.gz"), file("${baseName}.unmapped_R2.fastq.gz") into ch_surecall

    module 'java/openjdk-1.14.02' 
    module 'samtools'
    """
    bash /projects/vh83/local_software/agent3.0/agent.sh trim -fq1 ${R1} -fq2 ${R2} -v2 -out ./${baseName}.unmapped
    
    """
    
}

process alignBwa {
    
    label 'bwa'
    
    input:
        set baseName, file(R1), file(R2) from ch_surecall
    output:
        set baseName, file("${baseName}.aligned.bam") into ch_pipedBams, ch_forMetrics1

    publishDir path: './output/bams', mode: 'copy'

    module      bwaModule
    module      samtoolsModule
    
    
    script:
    """
    bwa mem -C -t ${task.cpus} $ref $R1 $R2 | samtools view -b - > ${baseName}.aligned.bam

    """
}


process setMateInfo {
    label 'genomics_3'

    input:
        set baseName, file(bam) from ch_pipedBams
    output:
        set baseName, file("${baseName}.aligned.matefixed.bam") into ch_mateFixed

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar SetMateInformation \
          -i $bam -r $ref -o ${baseName}.aligned.matefixed.bam
    """
}


process groupreadsByUmi {
    
    label 'genomics_3'

    input:
        set baseName, file(bam) from ch_mateFixed
    output:
        set baseName, file("${baseName}.piped.grouped.histogram.tsv"), file("${baseName}.piped.grouped.bam") into ch_umiGroupedBams
    
    publishDir path: './output/metrics/UMI/family_sizes', mode: 'copy', pattern: "*.tsv"

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar GroupReadsByUmi \
         -i ${bam} -f "${baseName}.piped.grouped.histogram.tsv" -o "${baseName}.piped.grouped.bam" -s Adjacency -e 1 
    """

}

process generateConsensusReads {
    
    label 'genomics_3'

    input:
        set baseName, file(hist), file(bam) from ch_umiGroupedBams
    output:
        set baseName, file("${baseName}.consensus.unmapped.bam") into ch_unmappedConsensusBams
    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar CallMolecularConsensusReads \
        --input $bam --output ${baseName}.consensus.unmapped.bam \
        --error-rate-post-umi 30 --min-reads 3
    """
}

process mapConsensusReads {
    
    label 'bwa'

    input:
        set baseName, file(bam) from ch_unmappedConsensusBams
    output:
        set baseName, file("${baseName}.consensus.aligned.bam") into ch_mappedConsensusBams, ch_forMetrics2
    publishDir path: './output/bams', mode: 'copy'

    module 	    bwaModule
    
    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR $tmp_dir | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin > ${baseName}.temp.bam
    samtools sort -o ${baseName}.consensus.aligned.bam ${baseName}.temp.bam
    """
}

process indexBam {
    
    label 'small_1'

    input:
        set baseName, file(bam) from ch_mappedConsensusBams
    output:
        set baseName, file(bam), file("${baseName}.consensus.aligned.bam.bai") into ch_indexedConsensusBams
    publishDir path: './output/bams', mode: 'copy'
 
    module      'samtools'
 
    script:
    """
    samtools index $bam ${baseName}.consensus.aligned.bam.bai
    """

}
