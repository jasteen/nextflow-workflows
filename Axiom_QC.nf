#!/usr/bin/env nextflow

chip_library_path = file("/projects/vh83/reference/axiom/r2/")
input_path = file("./cels/")

Channel
    .fromPath("${input_path}/*.cel")
    .into {ch_celList; ch_cels}


ch_celList.map { it -> it.name }
       .collectFile(name: 'cel.txt', newLine: true, seed: "cel_files")
       .set {ch_celList_QC}

//set the second to all the files

ch_cels
    .collect()
    .into {ch_cels_QC;ch_cels_GTQC;ch_cels_Summary;ch_cels_GT}


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

    label 'small_short'

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

/*
$CEL_LIST_INLIERS_2 is the path to a text file with the header “cel_files” and each subsequent 
row a path to each CEL file, which should all pass the DQC threshold AND the QC call_rate threshold
from the previous two steps.
*/
//generate summary callrates for all probesets
process summary_callrates {

  label 'small_3'

    input:
        file '*' from ch_GTQCout
        file '*' from ch_cels_Summary 

    output:
        set file("AxiomGT1.report.txt"), file("*") into ch_Summaryout
        
    publishDir path: './output/summary', mode: 'copy'

  script:
  """
  apt-genotype-axiom \
  --analysis-files-path ${chip_library_path} \
  --arg-file ${chip_library_path}/Axiom_ABC.r2.apt-genotype-axiom.AxiomCN_GT1.apt2.xml \
  --cel-files ./pass_GT_cel_list.txt \
  --out-dir ./ \
  --log-file ./apt2-axiom.log
  """
}

//need to make sure we capture the right things from the previous step.

process CNV {

  label 'small_3'

    input:
        set file(axoiom), file ('*') from ch_Summaryout
    output:
        file '*' into ch_CNVout
    
    publishDir path: './output/cn', mode: 'copy'

  script:
  """
  apt-copynumber-axiom-cnvmix \
  --analysis-files-path ${chip_library_path} \
  --arg-file ${chip_library_path}/Axiom_ABC.r2.apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml \
  --mapd-max 0.35 \
  --waviness-sd-max 0.1 \
  --summary-file ./AxiomGT1.summary.a5 \
  --report-file ./AxiomGT1.report.txt \
  --out-dir ./ \
  --log-file ./apt-copynumber-axiom.log
  """
}

process run_finalGT {
   
  label 'small_3'

    input:
        file '*' from ch_CNVout
        file '*' from ch_cels_GT

    output:
        file '*' into ch_GTout
    
    publishDir path: './output/genotypes', mode: 'copy'


  script:
  """
  apt-genotype-axiom \
  --copynumber-probeset-calls-file ./AxiomCNVMix.cnpscalls.txt \
  --analysis-files-path ${chip_library_path} \
  --arg-file ${chip_library_path}/Axiom_ABC_96orMore_Step2.r2.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
  --dual-channel-normalization true \
  --cel-files ./pass_GT_cel_list.txt \
  --out-dir ./ \
  --batch-folder ./ \
  --log-file ./apt2-axiom.log \
  --allele-summaries true \
  --write-models
  """
}

/*
process run_snpQC {
  script:
  """
  ps-metrics \
  --posterior-file $OUTDIR/genotypes/AxiomGT1.snp-posteriors.txt \
  --multi-posterior-file $OUTDIR/genotypes/AxiomGT1.snp-posteriors.multi.txt \
  --batch-folder $OUTDIR/genotypes \
  --summary-file $OUTDIR/genotypes/AxiomGT1.summary.txt \
  --report-file $OUTDIR/genotypes/AxiomGT1.report.txt \
  --special-snps $AXIOM_LIB_PATH/Axiom_PMDA.r6.specialSNPs \
  --use-multi-allele true \
  --y-restrict 0.2 \
  --metrics-file $OUTDIR/SNPolisher/metrics.txt \
  --multi-metrics-file $OUTDIR/SNPolisher/metrics.multi.txt \
  --log-file $OUTDIR/SNPolisher/ps_metrics.log
  """
}

process run_snpClassification {
  script:
  """
   ps-classification \
  --metrics-file $OUTDIR/SNPolisher/metrics.txt \
  --multi-metrics-file $OUTDIR/SNPolisher/metrics.multi.txt \
  --ps2snp-file $AXIOM_LIB_PATH/Axiom_PMDA.r6.ps2snp_map.ps \
  --species-type Human \
  --cr-cutoff 95 \
  --fld-cutoff 3.6 \
  --het-so-cutoff -0.1 \
  --het-so-XChr-cutoff -0.1 \
  --het-so-otv-cutoff -0.3 \
  --hom-ro-1-cutoff 0.6 \
  --hom-ro-2-cutoff 0.3 \
  --hom-ro-3-cutoff -0.9 \
  --hom-ro true \
  --hom-het true \
  --num-minor-allele-cutoff 2 \
  --hom-ro-hap-1-cutoff 0.1 \
  --hom-ro-hap-1-XChr-cutoff 0.1 \
  --hom-ro-hap-1-MTChr-cutoff 0.4 \
  --hom-ro-hap-2-cutoff -0.9 \
  --hom-ro-hap-2-XChr-cutoff 0.05 \
  --hom-ro-hap-2-MTChr-cutoff 0.2 \
  --hom-hap-X-cutoff -1 \
  --hom-hap-Y-lower-cutoff 1 \
  --hom-hap-Y-upper-cutoff 1 \
  --CN0-hap-X-cutoff -1 \
  --CN0-hap-Y-cutoff -1 \
  --CN0-dip-X-cutoff -1 \
  --CN0-dip-Y-cutoff -1 \
  --aaf-XChr-cut 0.36 \
  --fld-XChr-cut 4 \
  --homfld-XChr-cut 6.5 \
  --homfld-YChr-cut 6.5 \
  --min-YChr-samples-cut 5 \
  --sign-diff-hom-1-cutoff 0.4 \
  --sign-diff-hom-2-cutoff 0.4 \
  --min-mean-cp2-cutoff 8.5 \
  --max-mean-cp2-cutoff 15 \
  --priority-order PolyHighResolution,NoMinorHom,OTV,MonoHighResolution,CallRateBelowThreshold \
  --recommended PolyHighResolution,NoMinorHom,MonoHighResolution,Hemizygous \
  --use-multi-allele true \
  --output-dir $OUTDIR/SNPolisher \
  --log-file $OUTDIR/SNPolisher/ps_classification.log
  """
}
*/
