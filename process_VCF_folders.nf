#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
panel_int      = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.interval_list')
padded_int     = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.interval_list')
panel_bed      = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.bed')
padded_bed     = file('/projects/vh83/reference/genomes/b37/accessory_files/intervals_Broad.human.exome.b37.padded.bed')
tmp_dir        = file('/scratch/vh83/tmp/')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("/projects/vh83/reference/genomes/b37/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
knownIndels      = file("/projects/vh83/reference/genomes/b37/accessory_files/1000G_phase1.indels.b37.vcf")
dbSNP            = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")
ceu_mergeGvcf    = file("/projects/vh83/reference/genomes/b37/accessory_files/ugp-1k-backgrounds/CEU_mergeGvcf.vcf")
fin_mergeGvcf    = file("/projects/vh83/reference/genomes/b37/accessory_files/ugp-1k-backgrounds/FIN_mergeGvcf.vcf")
gbr_mergeGvcf    = file("/projects/vh83/reference/genomes/b37/accessory_files/ugp-1k-backgrounds/GBR_mergeGvcf.vcf")


header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

mills_grch37          = file("/projects/vh83/reference/genomes/b37/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
one_k_g_grch37_indels = file("/projects/vh83/reference/genomes/b37/accessory_files/1000G_phase1.indels.b37.vcf")
one_k_g_snps          = file("/projects/vh83/reference/genomes/b37/accessory_files/1000G_omni2.5.b37.vcf")
one_k_g_highconf_snps = file("/projects/vh83/reference/genomes/b37/accessory_files/1000G_phase1.snps.high_confidence.b37.vcf")
one_k_g_indels        = file("/projects/vh83/reference/genomes/b37/accessory_files/1000G_phase1.indels.b37.vcf")
hapmap                = file("/projects/vh83/reference/genomes/b37/accessory_files/hapmap_3.3.b37.vcf")
interval_grch37       = file("/projects/vh83/reference/genomes/b37/accessory_files/Broad.human.exome.b37.interval_list")
dbsnp_grch37          = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")

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



// Tools
picardJar          = '~/picard.jar'
gatkJar            = '/usr/local/gatk/3.7/bin/GenomeAnalysisTK.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
gatkModule         = 'gatk/3.7'
bcftoolsModule     = 'bcftools/1.8'
rModule            = 'R/3.5.1'          
fgbioJar           = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'


ch_gVcfs = Channel.fromPath("./variants/GATK/gvcf/*.vcf")
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
        file 'genotyped.vcf' into ch_genotypedGVCFsnp, ch_genotypedGVCFindel, ch_forCombining
    
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
                    -input $vcf --recal_file "output.recal_snp" --tranches_file "output.tranches_snp" \
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
        file(indel_recal) from ch_indelRecal
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

    publishDir path: './processed/gatk', mode: 'copy', pattern: "*.vcf"
    
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
        set sample, file("${sample}.vardict.sorted.vcf.gz") into ch_reheaderVCF

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
        set file(vcf), file("${sample}.vardict.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './output/vcf/UMI', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
    """
    bcftools index -f --tbi ${vcf} -o "${sample}.vardict.sorted.vcf.gz.tbi"
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

    publishDir path: './processed/vardict', mode: 'copy', pattern: "*.vcf"
    

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
                      -o "vardict.merged.vt.vep.vcf"
    """
}


