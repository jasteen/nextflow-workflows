#!/usr/bin/env nextflow

chip_library_path = file("/projects/vh83/reference/axiom/r2/")
input_path = file("./cels")


Channel
    .fromPath("${input_path}/*")
    .into {ch_celList; ch_cels}

ch_celList.map { it -> it.name }
       .collectFile(name: 'cel.txt', newLine: true, seed: "cel_files")
       .set {ch_celList_QC}

//set the second to all the files

ch_cels
    .collect()
    .into {ch_cels_QC;ch_cels_GTQC}


//run QC and remove samples with < 0.82 DQC
process runDQC {

    label 'small_short'

    input:
        file list from ch_celList_QC
        file '*' from ch_cels_QC   
    output:
        set file("raw_dqc.txt"), file("pass_DQC.txt") into ch_DCQCout

    publishDir path: './output/DQC', mode: 'copy'
    
    script:

    """
    apt-geno-qc --analysis-files-path ${chip_library_path} \
    --xml-file ${chip_library_path}/Axiom_ABC.r2.apt-geno-qc.AxiomQC1.xml \
    --cel-files cel.txt \
    --out-file raw_dqc.txt 
    awk 'BEGIN{FS=OFS="\t"}{if(NR == 1){print}else if(\$18 >= 0.82){ print \$1}}' raw_dqc.txt > pass_DQC.txt
    """
}

//run first pass genotype QC and remove samples with call rate <0.97
process runGTQC {

    label 'small_3'

    input:
        set file(raw_dqc), file (pass_dqc) from ch_DCQCout
        file '*' from ch_cels_GTQC   
    output:
        file '*' into ch_GTQCout
        
    publishDir path: './output/GTQC', mode: 'copy'
    
    script:
    """
    apt-genotype-axiom --analysis-files-path ${chip_library_path} \
    --arg-file  ${chip_library_path}/Axiom_ABC_96orMore_Step1.r2.apt-genotype-axiom.AxiomGT1.apt2.xml \
    --dual-channel-normalization true \
    --cel-files ./pass_DQC.txt \
    --table-output false \
    --out-dir . \
    --log-file ./apt2-axiom.log
    awk 'BEGIN{FS=OFS="\t"}{if(\$0 ~ /^#/){next}else if(\$1 ~ /cel_files/){print}else if(\$3 >= 97.0){ print \$1}}' AxiomGT1.report.txt > pass_GT.txt
    cut -f1 pass_GT.txt > pass_GT_cel_list.txt
    """
  }
  