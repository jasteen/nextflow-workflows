process merge_gvcf{
    """
    g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
        gatk_args = "-T CombineGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "{g_vcf_files} -o {vcf_out}".f
    """


}
process joint genotype {
    """
    gatk_args = "-T GenotypeGVCFs -R {reference} " \
                       "--disable_auto_index_creation_and_locking_when_reading_rods " \
                                          "--variant {merged_vcf} --out {vcf_out} " \
                  "--variant {CEU_mergeGvcf} --variant {GBR_mergeGvcf} "
    """

}
process snp_recal {

    """
    gatk_args = "-T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --num_threads {cores} " \
                    "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} " \
                    "-resource:omni,known=false,training=true,truth=true,prior=12.0 {one_k_g_snps} " \
                    "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {one_k_g_highconf_snps} " \
                    "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff " \
                    "-input {genotype_vcf} --recal_file {recal_snp} --tranches_file {tranches_snp} 
    """

}
process indel_recal {

    """
    java -Xmx6g -jar gatkJar -T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --num_threads {cores} " \
                    "-resource:mills,known=false,training=true,truth=true,prior=12.0 {mills_grch37} " \
                    "-resource:1000G,known=false,training=true,truth=true,prior=10.0 {one_k_g_indels} " \
                    "--maxGaussians 4 -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff " \
                    "-input {genotype_vcf} -recalFile {recal_indel} " \
                    "-tranchesFile {tranches_indel} -rscriptFile {indel_plots} " \
                    " -mode INDEL"
    """
}

process apply_indel{
    """
    gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
#                    "-R {reference} --ts_filter_level 99.0 --excludeFiltered --num_threads {cores} " \
#                    "-input {genotype_vcf} -recalFile {recal_indel} -tranchesFile {tranches_indel} " \
#                    "-mode INDEL -o {vcf_out}"
    """


}
process apply_snp {

 """
 gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
#                    "-R {reference} --ts_filter_level 99.5 --excludeFiltered --num_threads {cores} " \
#                    "-input {genotype_vcf} -recalFile {recal_snp} -tranchesFile {tranches_snp} " \
#                    "-mode SNP -o {vcf_out}"
"""
}


--interval_padding 50 
-L ${panel_int}
