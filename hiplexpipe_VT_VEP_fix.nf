#!/usr/bin/env nextflow

// Required Inputs
params.input_file = ""
params.reference = ""
//set up for multiple regerence possibilities.

if(params.reference == "hg19"){
    //HG19 reference for aspree stuff
    refFolder      = file("/projects/vh83/reference/genomes/hg19")
    refBase        = "$refFolder/ucsc.hg19"
    ref            = file("${refBase}.fasta")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")

}else if(params.reference == "hg38"){
    refFolder      = file("/projects/vh83/reference/genomes/hg38/heng_li_recomended")
    refBase        = "$refFolder/GCA_000001405.15_GRCh38_no_alt_analysis_set"
    ref            = file("${refBase}.fna")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")

}else{
    
    refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index")
    refBase          = "$refFolder/human_g1k_v37_decoy"
    ref              = file("${refBase}.fasta")
    refDict          = file("${refBase}.dict")
    refFai           = file("${refBase}.fasta.fai")
}


//Annotation resources
dbsnp_b37       = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")
other_vep       = file("/usr/local/vep/90/ensembl-vep/cache")


//extra stuff, may not be needed in this script
tmp_dir        = file('/scratch/vh83/tmp/')

// Creating channel from input directory
ch_VCF = Channel.fromPath(params.input_file).map {file -> tuple(file.baseName, file)} 

process vt_decompose_normalise {

    label 'genomics_1'
        
    input:
        set baseName, file(vcf) from ch_VCF
    output:
        set baseName, file("${baseName}.vt.vcf.gz") into ch_vtDecomposeVCF

    publishDir path: './', mode: 'copy'

    module      'vt'

    script:
    """
    vt decompose -s $vcf | vt normalize -n -r $ref -o ${baseName}.vt.vcf.gz -
    """
}

process apply_vep {

    label 'vep'

    input:
        set baseName, file(vcf) from ch_vtDecomposeVCF
    output:
        file("${baseName}.vt.vep.vcf") into ch_vepVCF

    publishDir path: './', mode: 'copy'

    module      'vep/90'

    script:
    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o ${baseName}.vt.vep.vcf
    """
}
