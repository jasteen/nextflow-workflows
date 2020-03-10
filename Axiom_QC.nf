#!/usr/bin/env nextflow

chip_library_path = file("/projects/vh83/reference/axiom")
input_path = file("./cels/")

Channel
    .fromPath("${input_path}/*.cel")
    .into {ch_celList; ch_cels}


ch_celList.map { it -> it.name }
       .collectFile(name: 'cel.txt', newLine: true, seed: "cel_files")
       .into {ch_celList_QC; ch_celList_GT}

//set the second to all the files

ch_cels
    .collect()
    .into {ch_cels_QC; ch_cels_GT}

process runQC {

    label 'small_short'

    input:
        file list from ch_celList_QC
        file '*' from ch_cels_QC   
    output:
        file("qc.txt") into ch_output

    publishDir path: './output', mode: 'copy'
    
    script:

    """
    apt-geno-qc --analysis-files-path ${chip_library_path} \
    --xml-file ${chip_library_path}/Axiom_ABC.r1.apt-geno-qc.AxiomQC1.xml \
    --cel-files cel.txt \
    --out-file qc.txt 
    """
}

process runGT {

    label 'medium_6h'

    input:
        file list from ch_celList_GT
        file '*' from ch_cels_GT   
    
    publishDir path: './output', mode: 'copy'
    
    script:
    """
    apt-probeset-genotype \
    --analysis-files-path ${chip_library_path} \
    --xml-file ${chip_library_path}/Axiom_ABC_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
    --out-dir ./out \
    --cel-files cel.txt
    """
  }



