#!/usr/bin/env nextflow

// Required Inputs
params.input_file = ""
params.reference == ""
//set up for multiple regerence possibilities.

if(params.reference == "hg19"){
    //HG19 reference for aspree stuff
    refFolder      = file("/projects/vh83/reference/genomes/hg19")
    refBase        = "$refFolder/ucsc.hg19"
    ref            = file("${refBase}.fasta")
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
vep_brcaex      = file("/projects/vh83/reference/annotation_databases/BRCA-Exchange/BRCA-exchange_accessed-180118/BRCA-exchange_accessed-180118.sort.vcf.gz")
vep_gnomad      = file("/projects/vh83/reference/annotation_databases/gnomAD/gnomad.exomes.r2.0.2.sites/gnomad.exomes.r2.0.2.sites.vcf.gz")
vep_revel       = file("/projects/vh83/reference/annotation_databases/REVEL/REVEL-030616/revel_all_chromosomes.vcf.gz")
vep_maxentscan  = file("/projects/vh83/reference/annotation_databases/MaxEntScan/MaxEntScan_accessed-240118")
vep_exac        = file("/projects/vh83/reference/annotation_databases/ExAC/ExAC_nonTCGA.r0.3.1/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz")
vep_dbnsfp      = file("/projects/vh83/reference/annotation_databases/dbNSFP/dbNSFPv2.9.3-VEP/dbNSFP-2.9.3.gz")
vep_dbscsnv     = file("/projects/vh83/reference/annotation_databases/dbscSNV/dbscSNV1.0-VEP/dbscSNV.txt.gz")
vep_cadd        = file("/projects/vh83/reference/annotation_databases/CADD/CADD-v1.3/1000G_phase3.tsv.gz")

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
                      --custom $vep_brcaex,brcaex,vcf,exact,0,Clinical_significance_ENIGMA,Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
                      --custom $vep_gnomad,gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
                      --custom $vep_revel,RVL,vcf,exact,0,REVEL_SCORE \
                      --plugin MaxEntScan,$vep_maxentscan \
                      --plugin ExAC,$vep_exac,AC,AN \
                      --plugin dbNSFP,$vep_dbnsfp,REVEL_score,REVEL_rankscore \
                      --plugin dbscSNV,$vep_dbscsnv \
                      --plugin CADD,$vep_cadd \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o ${baseName}.vt.vep.vcf
    """
}
