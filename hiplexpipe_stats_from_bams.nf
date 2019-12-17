#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files


intervalFile     = file("/projects/vh83/reference/hiplex_complexo/complexo_nthl1.bed")

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

// Tools
picardJar      = '~/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'
bamclipper_exe = '/projects/vh83/local_software/bamclipper/bamclipper.sh'

// Creating channel from input directory
ch_inputFiles1 = Channel.fromFilePairs("./bams/genesis/*.sort.hq.{bam,bam.bai}", flat: true)
ch_inputFiles2 = Channel.fromFilePairs("./bams/HVIO/*.hq.sorted.{bam,bam.bai}", flat: true)
ch_inputFiles3 = Channel.fromFilePairs("./bams/ABCFS/*.sort.hq.{bam,bam.bai}", flat: true)
ch_inputFiles4 = Channel.fromFilePairs("./bams/ABCFS/halo/*.sorted.locatit.{bam,bam.bai}", flat: true)

ch_inputFiles5 = ch_inputFiles1.join(ch_inputFiles2)
ch_inputFiles6 = ch_inputFiles5.join(ch_inputFiles3)
ch_inputFiles7 = ch_inputFiles6.join(ch_inputFiles4)

ch_inputFiles7.into{ch_mappedBam1;ch_mappedBam2;ch_mappedBam3;ch_mappedBam4;ch_mappedBam5}

/*
Stats Generation Section
*/

process IntersectBed {

    label 'medium_1h'

    input:
        set sample, file(bam), file(bai) from ch_mappedBam2
    output:
        set sample, file("${sample}.intersectbed.bam") into ch_intersectBam
    
    script:
    """
    module load bedtools/2.27.1-gcc5
    intersectBed -abam ${bam} -b ${intervalFile} > ${sample}.intersectbed.bam
    """
}

process CoverageBed {

    label 'medium_1h'

    input:
        set sample, file(bam), file(bai) from ch_mappedBam3
    output:
        set sample, file("${sample}.bedtools_hist_all.txt"), file("${sample}.bedtools_amplicon.txt") into ch_bedtools
    
    script:
    """
    module load bedtools/2.27.1-gcc5
    coverageBed -b ${bam} -a ${intervalFile} \
        -sorted -hist -g ${genome_file} | \
        grep all > "${sample}.bedtools_hist_all.txt"
    bedtools coverage -f 5E-1 -a ${intervalFile} -b ${bam} | sed "s/\$/    ${sample}/g" > "${sample}.bedtools_amplicon.txt"
    """
}

process ReadsMapped {
    
    label 'medium_1h'

    input:
        set sample, file(bam), file(bai) from ch_mappedBam4
    output:
        set sample, file("${sample}.mapped_to_genome.txt") into ch_onGenome

    module      'samtools/1.9'

    script:
    """
    samtools view -c -F4 ${bam} > "${sample}.mapped_to_genome.txt"
    """
}

process ReadsTotal {

    label 'medium_1h'

    input:
        set sample, file(bam), file(bai) from ch_mappedBam5
    output:
        set sample, file("${sample}.total_raw_reads.txt") into ch_onTotal

    module      'samtools/1.9'
   
    script:
    """
    samtools view -c ${bam} > "${sample}.total_raw_reads.txt"
    """
}
    
process TargetMapped {
    
    label 'medium_1h'

    input:
        set sample, file(bam) from ch_intersectBam
    output:
        set sample, file("${sample}.mapped_to_target.txt") into ch_onTarget

    module      'samtools/1.9'

    script:
    """
    samtools view -c -F4 ${bam} > ${sample}.mapped_to_target.txt
    """
}

ch_final = ch_bedtools.join(ch_onGenome)
ch_final2 = ch_final.join(ch_onTarget)
ch_final3 = ch_final2.join(ch_onTotal)

process collateData {
    
    label 'medium_1h'

    input:
        set sample, file(bedtools), file(bedtools2), file(onGenome), file(onTarget), file(onTotal) from ch_final3
    output:
        set sample, file("${sample}_summary_coverage.txt") into ch_out

    script:
    """
    module purge
    module load R/3.5.1
    Rscript --vanilla /projects/vh83/pipelines/code/modified_summary_stat.R \
            ${bedtools} \
            ${onGenome} \
            ${onTarget} \
            ${onTotal} \
            ${sample} \
            "${sample}_summary_coverage.txt"
    """
}

ch_out.map{a,b -> b}.collect().set{ch_out2}

process catStats {

    label 'small_1'

    input:
        file(stats) from ch_out2
    output:
        file("project_summary.txt") into ch_out3
    
    publishDir path: './metrics/', mode: 'copy'

    script:
    """
    cat ${rheader} ${stats} > "project_summary.txt"
    """

}


