#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bowtie2_2.2.9_index")
inputDirectory = file('./fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files

vardictBed       = file("/projects/vh83/reference/CONFIRM/vardict/CONFIRM_8Col.bed")
intervalFile     = file("/projects/vh83/reference/confirm/CONFIRM20_PMS2_67810_nobreak.final.numsort_b37.bed")
restrictedBed    = file("/projects/vh83/reference/confirm/CONFIRM20_PMS2_67810_nobreak.final.numsort_b37.bed")
primer_bedpe_file= file("/projects/vh83/reference/confirm/CONFIRM20_PMS2_67810_nobreak.final.bamclipper.bedpe.b37.tsv")

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
genome_file      = file("/projects/vh83/reference/genomes/b37/accessory_files/human_g1k_v37_decoy_GenomeFile.txt")

// Creating channel from input directory
ch_inputFiles = Channel.fromFilePairs("${inputDirectory}/*_R{1,2}_001.fastq.gz")


process align_bwa {

    label 'bwa_small'

    input:
        set baseName, file(fastqs) from ch_inputFiles
    output:
        set baseName, file("${baseName}.bowtie.sorted.bam"), file("${baseName}.bowtie.sorted.bam.bai") into ch_mappedBams

    publishDir path: './bam_out', mode: 'copy'

    module      'bowtie2/2.2.9'
    module      'samtools/1.9'

    """
    bowtie2 -p ${task.cpus} -x $ref --trim5 0 --trim3 60 -1 ${fastqs[0]} -2 ${fastqs[1]} \
        | samtools view -bh - \
        | samtools sort -@ ${task.cpus} -o "${baseName}.bowtie.sorted.bam"
    samtools index "${baseName}.bowtie.sorted.bam" "${baseName}.bowtie.sorted.bam.bai"
    """
}