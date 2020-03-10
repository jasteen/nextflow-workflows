#!/usr/bin/env nextflow

chip_library_path = file("/projects/vh83/reference/axiom/")
input_path = file("./cels/")

Channel
    .fromPath("${input_path}/*.cel")
    .into {ch_celList; ch_cels}


ch_celList.map { it -> it[1].name }
       .collectFile(name: 'cel.txt', newLine: true, seed: "cel_files")
       .set {ch_celList_f}

//set the second to all the files

ch_cels
    .collect()
    .set {ch_cels}

process runQC {

    label 'medium_6h'

    input:
        file list from ch_celList_f
        file '*' from ch_cels   

    publishDir path: './output', mode: 'copy'
    
    script:

    """
    apt-geno-qc --analysis-file-path ${chip_library_path} \
    --xml-file Axiom_ABC.r1.apt-geno-qc.AxiomQC1.xml \
    --cel-files cel.txt
    --out qc.txt 
    """
}

