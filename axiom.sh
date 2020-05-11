!#/bin/bash

ls ./cels/*.cel > cel_list.txt
ORIGINAL_CEL_FILES=cel_list.txt
AXOIOM_LIB_PATH=/projects/vh83/reference/axiom/r2/
mkdir ./bash
OUTDIR=./bash

apt-geno-qc \
	--analysis-files-path $AXIOM_LIB_PATH \
	--xml-file $AXIOM_LIB_PATH/Axiom_ABC.r2.apt-geno-qc.AxiomQC1.xml \
	--cel-files $ORIGINAL_CEL_FILES \
	--out-dir $OUTDIR/qc \
	--out-file $OUTDIR/qc/apt-geno-qc.report.txt \
	--log-file $OUTDIR/qc/apt-geno-qc.log

awk 'BEGIN{FS=OFS="\t"}{if(NR == 1){print}else if(\$18 >= 0.82){ print \$1}}' \
         ./$OUTPUT/apt-geno-qc.report.txt > ./$OUTPUT/pass_DQC.txt

CEL_LIST_INLIERS=./$OUTPUT/pass_DQC.txt

apt-genotype-axiom \
	--analysis-files-path $AXIOM_LIB_PATH \
	--arg-file $AXIOM_LIB_PATH/Axiom_ABC_96orMore_Step1.r2.apt-genotype-axiom.AxiomGT1.apt2.xml \
	--dual-channel-normalization true \
	--cel-files $CEL_LIST_INLIERS_1 \
	--table-output false \
	--out-dir $OUTDIR/qc \
	--log-file $OUTDIR/qc/apt2-axiom.log

awk 'BEGIN{FS=OFS="\t"}{if(\$0 ~ /^#/){next}else if(\$1 ~ /cel_files/){print}else if(\$3 >= 97.0){ print \$1}}' ./$OUTPUT/qc/AxiomGT1.report.txt > ./$OUTPUT/qc/pass_GT.txt
    cut -f1 ./$OUTPUT/qc/pass_GT.txt > ./$OUTPUT/qc/pass_GT_cel_list.txt

CEL_LIST_INLIERS_2=./$OUTPUT/qc/pass_GT_cel_list.txt

apt-genotype-axiom \
	--analysis-files-path $AXIOM_LIB_PATH \
	--arg-file $AXIOM_LIB_PATH/Axiom_ABC.r2.apt-genotype-axiom.AxiomCN_GT1.apt2.xml \
	--cel-files $CEL_LIST_INLIERS_2 \
	--out-dir $OUTDIR/summary \
	--log-file $OUTDIR/summary/apt2-axiom.log

apt-copynumber-axiom-cnvmix \
	--analysis-files-path $AXIOM_LIB_PATH \
	--arg-file $AXIOM_LIB_PATH/Axiom_ABC.r2.apt-copynumber-axiom-cnvmix.AxiomCNVmix.apt2.xml \
	--mapd-max 0.35 \
	--waviness-sd-max 0.1 \
	--summary-file $OUTDIR/summary/AxiomGT1.summary.a5 \
	--report-file $OUTDIR/summary/AxiomGT1.report.txt \
	--out-dir $OUTDIR/cn \
	--log-file $OUTDIR/cn/apt-copynumber-axiom.log

apt-genotype-axiom \
	--copynumber-probeset-calls-file $OUTDIR/cn/AxiomCNVMix.cnpscalls.txt \
	--analysis-files-path $AXIOM_LIB_PATH \
	--arg-file $AXIOM_LIB_PATH/Axiom_ABC_96orMore_Step2.r2.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
	--dual-channel-normalization true \
	--cel-files $CEL_LIST_INLIERS_2 \
	--out-dir $OUTDIR/genotypes \
	--batch-folder $OUTDIR/genotypes \
	--log-file $OUTDIR/genotypes/apt2-axiom.log \
	--allele-summaries true \
	--write-models



ps-metrics \
	--posterior-file $OUTDIR/genotypes/AxiomGT1.snp-posteriors.txt \
	--multi-posterior-file $OUTDIR/genotypes/AxiomGT1.snp-posteriors.multi.txt \
	--batch-folder $OUTDIR/genotypes \
	--summary-file $OUTDIR/genotypes/AxiomGT1.summary.txt \
	--report-file $OUTDIR/genotypes/AxiomGT1.report.txt \
	--special-snps $AXIOM_LIB_PATH/Axiom_ABC.r2.specialSNPs \
	--use-multi-allele true \
	--y-restrict 0.2 \
	--metrics-file $OUTDIR/SNPolisher/metrics.txt \
	--multi-metrics-file $OUTDIR/SNPolisher/metrics.multi.txt \
	--log-file $OUTDIR/SNPolisher/ps_metrics.log


ps-classification \
	--metrics-file $OUTDIR/SNPolisher/metrics.txt \
	--multi-metrics-file $OUTDIR/SNPolisher/metrics.multi.txt \
	--ps2snp-file $AXIOM_LIB_PATH/Axiom_ABC.r2.ps2snp_map.ps \
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

