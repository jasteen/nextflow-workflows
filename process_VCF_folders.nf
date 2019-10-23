

ch_gVcfs = Channel.fromPath("./variants/GATK/*.vcf")
ch_vardictVCFS = Channel.fromPath("./variants/vardict/*.vcf").map{file -> tuple(file.name.take(file.name.lastIndexOf('.')), file)}

//GATK section

process mergeGVCFS {
    
    label 'gatk_unknown'
    echo 'true'
    input:
        file(vcf) from ch_gVcfs.collect()
    output:
        file("combined.g.vcf") into ch_combinedGVCF

    module gatkModule

    script:

    myfiles = vcf.join(' -V ')
    
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T CombineGVCFs -R ${ref} \
                  --disable_auto_index_creation_and_locking_when_reading_rods \
                  -V $myfiles -o "combined.g.vcf"
    """
}

process genotypeGVCF {
    
    label 'gatk_unknown'
    
    input:
        file(vcf) from ch_combinedGVCF
    output:
        file("genotyped.vcf") into ch_genotypedGVCFsnp, ch_genotypedGVCFindel, ch_forCombining
    
    module gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T GenotypeGVCFs -R ${ref} \
                    --disable_auto_index_creation_and_locking_when_reading_rods \
                    --variant $vcf --out "genotyped.vcf" \
                    --variant ${ceu_mergeGvcf} --variant ${gbr_mergeGvcf} --variant ${fin_mergeGvcf}
    """
}

process snpRecalibrate {

    label 'gatk_unknown'
    
    input:
        file(vcf) from ch_genotypedGVCFsnp
    output:
        set file(vcf), file("output.recal_snp"), file("output.tranches_snp"), file("output.plots_snp") into ch_applysnpRecal

    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T VariantRecalibrator \
                    --disable_auto_index_creation_and_locking_when_reading_rods \
                    -R $ref --num_threads ${task.cpus} \
                    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
                    -resource:omni,known=false,training=true,truth=true,prior=12.0 $one_k_g_snps \
                    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $one_k_g_highconf_snps \
                    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
                    -input $vcf --recal_file "output.recal_snp" --tranches_file $output.tranches_snp \
                    -rscriptFile "output.plots_snp" -mode SNP
    """
}

process indelRecalibrate {
    
    label 'gatk_unknown'
    
    input:
    file(vcf) from ch_genotypedGVCFindel
    output:
    set file(vcf), file("output.recal_indel"), file("output.tranches_indel"), file("output.plots_indel") into ch_applyindelRecal
    
    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar  -T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods\
                    -R $ref --num_threads ${task.cpus} \
                    -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills_grch37 \
                    -resource:1000G,known=false,training=true,truth=true,prior=10.0 $one_k_g_indels \
                    --maxGaussians 4 -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
                    -input $vcf -recalFile "output.recal_indel" \
                    -tranchesFile "output.tranches_indel" -rscriptFile "output.plots_indel" \
                    -mode INDEL
    """
}

process applySNPrecal{

    label 'gatk_unknown'
    
    input:
        set file(vcf), file(recal), file(tranch), file(plots) from ch_applysnpRecal
    output:
        file("snp_recal.vcf") into ch_snpRecal
    
    module      gatkModule

    script:
    """   
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T ApplyRecalibration \
                    --disable_auto_index_creation_and_locking_when_reading_rods \
                    -R $ref --ts_filter_level 99.5 --excludeFiltered --num_threads ${task.cpus} \
                    -input $vcf -recalFile $recal -tranchesFile $tranches \
                    -mode SNP -o "snp_recal.vcf"
    """
}

process applyINDELrecal{
    
    label 'gatk_unknown'

    input:
        set file(vcf), file(recal), file(tranch), file(plots) from ch_applyindelRecal
    output:
        file("indel_recal.vcf") into ch_indelRecal    
    
    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods \
                    -R $ref --ts_filter_level 99.0 --excludeFiltered --num_threads ${task.cpus} \
                    -input $vcf -recalFile $recal -tranchesFile $tranches \
                    -mode INDEL -o "indel_recal.vcf"

    """
}

process combineAllRecal {
    
    label 'gatk_unknown'
    
    input:
        file(snp_recal) from ch_snpRecal
        file(indel_recal) from ch_yindelRecal
    output:
        file("recalibrated.vcf") into ch_finalGATKvcf
    //publishDir path: './output/GATK', mode: 'copy'

    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T CombineVariants \
                    -R $ref --disable_auto_index_creation_and_locking_when_reading_rods \
                    --num_threads ${task.cpus} --genotypemergeoption UNSORTED --variant $snp_recal \
                    --variant $indel_recal -o "recalibrated.vcf"
    """
}

process apply_vep_GATK {
    
    label 'vep'

    input:
        file(vcf) from ch_finalGATKvcf
    output:
        file("gatk.merged.vt.vep.vcf") into ch_annotateGATK

    publishDir path: './output/vcf/', mode: 'copy', pattern: "*.vcf"
    
    module      'vep/90'

    script:

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick \
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
                      -o "gatk.merged.vt.vep.vcf"
    """
}



//Vardict section

process reheaderUMIVCF {
    
    label 'small_1'

    input:
        set sample, file(vcf) from ch_vardictVCFS
    output:
        set sample, file("*.vcf.gz") into ch_reheaderVCF

    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    module      'bcftools/1.8'

    script:
  
    """
    bcftools annotate -h ~/vh83/reference/genomes/b37/vcf_contig_header_lines.txt -O v ${vcf} | \
        bcftools sort -o ${sample}.vardict.sorted.vcf.gz -O z -
    """
}

process indexVCFS {

    label 'small_1'

    input:
        set sample, file(vcf) from ch_reheaderVCF
    output:
        set sample, file(vcf), file("${sample}.vardict.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './output/vcf/UMI', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
  
    """
    bcftools index -f --tbi ${vcf} -o ${sample}.UMI.reheader.sorted.vcf.gz.tbi
    """
}

//duplicate ch_indexedVCF
ch_indexedVCF
    .into{ch_premergelist;ch_premerge_files}
//set one version to a list of filenames of the VCF
ch_premergelist.map { it -> it[0].name }
       .collectFile(name: 'list.txt', newLine: true)
       .set {ch_premerge_list_f}
//set the second to all the files
ch_premerge_files
    .collect()
    .set {ch_premerge_all_files}

//feed both to the merge so that the indexes are available to bcftools
process mergeVCFS {
    
    label 'small_1'
    
    publishDir './variants_merged/', mode: 'copy'
    
    input:
    file list from ch_premerge_list_f
    file '*' from ch_premerge_all_files
    
    output:
    file("vardict.merged.vcf.gz") into ch_mergedfinalVCF

    module     'bcftools/1.8'
    
    script: 
    
    """
    bcftools merge -m none -O z -o "vardict.merged.vcf.gz" -l list.txt
    """
}

process vt_decompose_normalise {
        
    label 'medium_6h'

    input:
        file(vcf) from ch_mergedfinalVCF
    output:
        file("vardict.merged.vt.vcf.gz") into ch_vtDecomposeVCF

    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    module      'vt'
    
    script:
    """
    vt decompose -s $vcf | vt normalize -r $ref -o "vardict.merged.vt.vcf.gz" -
    """
}

process apply_vep {
    
    label 'vep'

    input:
        file(vcf) from ch_vtDecomposeVCF
    output:
        file("vardict.merged.vt.vep.vcf") into ch_vepVCF

    publishDir path: './output/vcf/', mode: 'copy', pattern: "*.vcf"
    

    module      'vep/90'

    script:

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick \
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
                      -o "${sample}.UMI.reheader.sorted.vt.vep.vcf"
    """
}


