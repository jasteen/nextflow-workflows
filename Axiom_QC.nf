#!/usr/bin/env nextflow

chip_library_path = file("/projects/vh83/reference/axiom")
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
    .into {ch_cels_QC;ch_cels_GT;ch_cels_Summary}


//run QC and remove samples with < 0.82 DQC
process runQC {

    label 'small_short'

    input:
        file list from ch_celList_QC
        file '*' from ch_cels_QC   
    output:
        set file("raw_qc.txt"), file("pass_QC.txt") into ch_output

    publishDir path: './output/QC', mode: 'copy'
    
    script:

    """
    apt-geno-qc --analysis-files-path ${chip_library_path} \
    --xml-file ${chip_library_path}/Axiom_ABC.r1.apt-geno-qc.AxiomQC1.xml \
    --cel-files cel.txt \
    --out-file raw_qc.txt 
    awk 'BEGIN{FS=OFS="\t"}{if(NR == 1){print}else if(\$18 >= 0.82){ print \$1}}' raw_qc.txt > pass_QC.txt
    """
}

//run first pass genotype QC and remove samples with call rate <0.97
process runGTQC {

    label 'medium_6h'

    input:
        set file(raw_qc), file (pass_qc) from ch_output
        file '*' from ch_cels_GT   
    output:
        set file("AxiomGT1.report.txt"), file("apt2-axiom.log") into ch_GTQC_out
    
    publishDir path: './output/GTQC', mode: 'copy'
    
    script:
    """
    apt-genotype-axiom --analysis-files-path ${chip_library_path} \
  --arg-file  ${chip_library_path}/Axiom_ABC_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
  --dual-channel-normalization true \
  --cel-files ./pass_QC.txt \
  --table-output false \
  --out-dir . \
  --log-file ./apt2-axiom.log
    """
  }

/*
//generate summary callrates for all probesets
process summary_callrates {

  label 'medium_6h'

  script:
  """
  apt-genotype-axiom \
  --analysis-files-path $AXIOM_LIB_PATH \
  --arg-file $AXIOM_LIB_PATH/Axiom_PMDA.r6.apt-genotype-axiom.AxiomCN_PS1.apt2.xml \
  --cel-files $CEL_LIST_INLIERS_2 \
  --out-dir $OUTDIR/summary \
  --log-file $OUTDIR/summary/apt2-axiom.log
  """
}

  process CNV {
  script:
  """
  apt-copynumber-axiom-cnvmix \
  --analysis-files-path $AXIOM_LIB_PATH \
  --arg-file $AXIOM_LIB_PATH/Axiom_PMDA.r6.apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml \
  --mapd-max 0.35 \
  --waviness-sd-max 0.1 \
  --summary-file $OUTDIR/summary/AxiomGT1.summary.a5 \
  --report-file $OUTDIR/summary/AxiomGT1.report.txt \
  --out-dir $OUTDIR/cn \
  --log-file $OUTDIR/cn/apt-copynumber-axiom.log
  """
}

process run_finalGT {
  script:
  """
  apt-genotype-axiom \
  --copynumber-probeset-calls-file $OUTDIR/cn/AxiomCNVMix.cnpscalls.txt \
  --analysis-files-path $AXIOM_LIB_PATH \
  --arg-file $AXIOM_LIB_PATH/Axiom_PMDA_96orMore_Step2.r6.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
  --dual-channel-normalization true \
  --cel-files $CEL_LIST_INLIERS_2 \
  --out-dir $OUTDIR/genotypes \
  --batch-folder $OUTDIR/genotypes \
  --log-file $OUTDIR/genotypes/apt2-axiom.log \
  --allele-summaries true \
  --write-models
  """

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
