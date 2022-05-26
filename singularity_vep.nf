#!/usr/bin/env nextflow

// Required Inputs
params.reference = ""
params.input_file = ""


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


ch_VCF = Channel.fromPath(params.input_file).map {file -> tuple(file.baseName, file)}

process apply_vep {
    label 'vep_sing'

    input:
        set basename, file(vcf) from ch_mergedfinalVCF
    output:
        file("${basename}_vep.vcf") into ch_vepVCF

    publishDir path: './variants_merged', mode: 'copy'

    module 'singularity'

    script:
    """
    vep --cache --dir_cache $vep_cache \
                      --assembly GRCh38 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o ${basename}_vep.vcf
    """
}