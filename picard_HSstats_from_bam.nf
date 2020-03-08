input_path =file("./bams")
//panel_int      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Covered_UCSChg19.interval_list ')
//padded_int     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded_UCSChg19.interval_list')


panel_int      = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.interval_list')
padded_int     = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.interval_list')
panel_bed      = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.bed')
padded_bed     = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.padded.bed')
//Variables

refFolder      = file("/projects/vh83/reference/genomes/hg19/")
refBase          = "$refFolder/ucsc.hg19"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")



picardJar          = '~/picard.jar'


Channel
    .fromPath("${input_path}/*.bam")
    .map{ file -> tuple(file.name.take(file.name.lastIndexOf('.')), file) }
    .set { ch_1 }

process collectHSMetrics {
    label 'medium_6h'
    input:
        set sample, file(bam) from ch_1
    output:
        set sample, file("*.HSmetrics.txt"), file("*.perbase.txt"), file("*.pertarget.txt") into ch_metrics_unused2
    
    publishDir path: './output/metrics/coverage', mode: 'copy'

    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx6G -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${sample}.HSmetrics.txt" \
        -R ${ref} \
        -BI $panel_int \
        -TI $padded_int \
        --PER_BASE_COVERAGE "${sample}.perbase.txt" \
        --PER_TARGET_COVERAGE "${sample}.pertarget.txt"
    """
}


