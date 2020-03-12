input_path = "/scratch/uc23/hfettke/cfDNA_BAMs_41/extra_bams/postconnor_LOD"
bed_target=file("/scratch/uc23/jste0021/cellfree.bed")


//Variables
reference=file("/projects/uc23/reference/genomes/bwa_index/human_g1k_v37_decoy.fasta")
genome_file=file("/projects/vh83/reference/genomes/b37/accessory_files/human_g1k_v37_decoy_GenomeFile.txt")

// Global Resource Configuration Options
Channel
    .fromPath("${input_path}/*.bam")
    .map{ file -> tuple(file.name.take(file.name.lastIndexOf('.')), file) }
    .into { ch_1; ch_2; ch_3 }

process InstersectBed {
    label 'medium_1h'
    input:
        set sample, file(bam) from ch_1
    output:
        set sample, file("${sample}.intersectbed.bam") into ch_intersectBam
    


    script:
    """
    module load bedtools/2.27.1-gcc5
    intersectBed -abam ${bam} -b ${bed_target} > ${sample}.intersectbed.bam
    """
}

process CoverageBed {
    label 'medium_1h'
    input:
        set sample, file(bam) from ch_2
    output:
        set sample, file("${sample}.bedtools_hist_all.txt") into ch_bedtools

    script:
    """
    module load bedtools/2.27.1-gcc5
    coverageBed -b ${bam} -a ${bed_target} \
        -sorted -hist -g ${genome_file} | \
        grep all > "${sample}.bedtools_hist_all.txt"
    """
}
process ReadsMapped {
    label 'medium_1h'
    input:
        set sample, file(bam) from ch_3
    output:
        set sample, file("${sample}.mapped_to_genome.txt") into ch_onGenome


    module      'samtools/1.9'

    
    script:
    """
    samtools view -c -F4 ${bam} > "${sample}.mapped_to_genome.txt"
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

process collateData {
    label 'medium_1h'
    input:
        set sample, file(bedtools), file(onGenome), file(onTarget) from ch_final2
    output:
        set sample, file("${sample}_summary_coverage.txt") into ch_out
    
    publishDir path: './output/', mode: 'copy'

    script:
    """
    module purge
    module load R/3.5.1
    Rscript --vanilla /scratch/uc23/jste0021/stats.R \
            ${bedtools} \
            ${onGenome} \
            ${onTarget} \
            ${sample} \
            "${sample}_summary_coverage.txt"
    """
}

