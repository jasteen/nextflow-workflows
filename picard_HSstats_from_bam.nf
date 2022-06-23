inputDirectory = file('./bams')
params.reference = ""

params.panel_bed       = ""
params.padded_bed      = ""
params.panel_int       = ""
params.padded_int      = ""

//set up for multiple regerence possibilities.

if(params.reference == "hg19"){
    //HG19 reference for aspree stuff
    refFolder      = file("/projects/vh83/reference/genomes/hg19")
    refBase        = "$refFolder/ucsc.hg19"
    ref            = file("${refBase}.fasta")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")

}else if(params.reference == "hg38"){
    refFolder      = file("/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0")
    refBase        = "$refFolder/Homo_sapiens_assembly38"
    ref            = file("${refBase}.fasta")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")
    genome_file    = file("$refFolder/hg38.chrom.sizes")
    header         = file("$refFolder/hg38_vcf_header.txt")
    vep_cache      = file("/projects/vh83/reference/VEP_CACHE")


}else{
    
    refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index")
    refBase          = "$refFolder/human_g1k_v37_decoy"
    ref              = file("${refBase}.fasta")
    refDict          = file("${refBase}.dict")
    refFai           = file("${refBase}.fasta.fai")
}


picardJar      = '/usr/local/picard/2.19.0/bin/picard.jar'


Channel
    .fromPath("${inputDirectory}/*.bam")
    .map{ file -> tuple(file.name.take(file.name.lastIndexOf('.')), file) }
    .set { ch_1 }

process collectHSMetrics {

    label 'start_1_8_2h'

    input:
        set sample, file(bam) from ch_1
    output:
        set sample, file("*.HSmetrics.txt"), file("*.pertarget.txt") into ch_metrics
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    script:

    """
    module purge
    module load R/3.5.1
    module picard/2.19.0

    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${bam.baseName}.HSmetrics.txt" \
        -R ${ref} \
        -BI ${params.panel_int} \
        -TI ${params.padded_int} \
        --PER_TARGET_COVERAGE "${bam.baseName}.pertarget.txt"
    """
}


