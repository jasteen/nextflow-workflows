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
ceu_mergeGvcf    = file("reference/CEU_mergeGvcf.vcf")
fin_mergeGvcf    = file("reference/FIN_mergeGvcf.vcf")
gbr_mergeGvcf    = file("reference/GBR_mergeGvcf.vcf")

mills_grch37 = file("Mills_and_1000G_gold_standard.indels.b37.vcf")
one_k_g_grch37_indels = file("1000G_phase1.indels.b37.vcf")
one_k_g_snps = file("1000G_omni2.5.b37.vcf")
one_k_g_highconf_snps = file("1000G_phase1.snps.high_confidence.b37.vcf")
one_k_g_indels = file("1000G_phase1.indels.b37.vcf")
hapmap = file("hapmap_3.3.b37.vcf")
interval_grch37 = file("Broad.human.exome.b37.interval_list")
dbsnp_grch37 = file("dbsnp_138.b37.vcf")

// Tools
picardJar          = '~/picard.jar'
gatkJar            = '/usr/local/gatk/3.7/bin/GenomeAnalysisTK.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
gatkModule         = 'gatk/3.7'
bcftoolsModule     = 'bcftools/1.8'
rModule            = 'R/3.5.1'          
fgbioJar           = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'


//initial channel generation from fastq folder
ch_inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz", flat: true)

process createUnmappedBam {
    
    label '3h_6g' 

    //publishDir path: './output/intermediate', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2) from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped.bam") into ch_unmappedBams

    publishDir path: './output/intermediate', mode: 'copy'
    
    module      'fgbio'
    module      'java'
    
    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -XX:+AggressiveOpts -XX:+AggressiveHeap \
        -jar $fgbioJar FastqToBam --input $R1 $R2 --output "${baseName}.unmapped.bam" --read-structures +T +T \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}

process markAdaptors {

    label '3h_6g' 

    //publishDir path: './output/intermediate', mode: 'copy'

    input:
        set baseName, file(bam) from ch_unmappedBams
    output:
        set baseName, file("${baseName}.unmapped.marked.bam"),
                      file("${baseName}.unmapped.marked_metrics.tsv") into ch_markedBams
   
    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar MarkIlluminaAdapters \
        -INPUT $bam \
        -OUTPUT "${baseName}.unmapped.marked.bam" \
        -METRICS "${baseName}.unmapped.marked_metrics.tsv"
    """
}

process alignBwa {
    
    label 'bwa'

    input:
        set baseName, file(bam), file(metrics) from ch_markedBams
    output:
        set baseName, file("${baseName}.mapped.bam") into ch_mappedBams

    //publishDir path: './output/bams', mode: 'copy'

    module	    samtoolsModule
    module      bwaModule

    script:
    """
    set -o pipefail
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ '/dev/stdout' -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 \
        -INTERLEAVE true -NON_PF true -TMP_DIR "$tmp_dir" | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.mapped.bam" -R "$ref" -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -SO queryname -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process markDuplicatesPicard {
    
    label '3h_6g'

    input:
        set baseName, bam from ch_mappedBams 
    output:
        set baseName, file("${baseName}.mapped.marked.bam") into ch_markedBamFiles
        set baseName, file("${baseName}.markduplicates.metrics") into ch_metrics

    publishDir path: './output/metrics/markduplicates', mode: 'copy'

    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar MarkDuplicates \
        -INPUT $bam \
        -OUTPUT ${baseName}.mapped.marked.bam \
        -METRICS_FILE ${baseName}.markduplicates.metrics \
        -VALIDATION_STRINGENCY SILENT \
        -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        -ASSUME_SORT_ORDER queryname
    """
}

process sortBam {

    label '3h_6g'

    input:
        set baseName, file(markedBam) from ch_markedBamFiles
    output:
        set baseName, file("${baseName}.mapped.marked.sorted.bam") into ch_sortedBamFiles

    script:
    """
    java -Djava.io.tmpdir=$tmp_dir -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar SortSam \
        -INPUT $markedBam \
        -OUTPUT ${baseName}.mapped.marked.sorted.bam \
        -SORT_ORDER coordinate \
        -CREATE_INDEX false \
        -CREATE_MD5_FILE false \
        -MAX_RECORDS_IN_RAM 300000
    """
}

process indexBam {

    label '3h_6g'

    input:
        set baseName, file(bam) from ch_sortedBamFiles
    output:
        set baseName, file(bam), file("${baseName}.mapped.marked.sorted.bam.bai") into ch_forVARDICT, ch_forGATK, ch_forHSMetrics, ch_forMultipleMetrics, ch_forperBase
    
    publishDir path: './output/bams', mode: 'copy'

    module      samtoolsModule

    script:
    """
    samtools index $bam ${baseName}.mapped.marked.sorted.bam.bai
    """
}

ch_forperBase.into{ch_bamList;ch_bams}
//set one version to a list of filenames of the VCF
ch_bamList.map { it -> it[1].name }
       .collectFile(name: 'list.txt', newLine: true)
       .set {ch_bamList_f}
//set the second to all the files
ch_bams
    .collect()
    .set {ch_all_bams}

process generatePerbaseMetrics {
    
    label 'medium_6h'

    input:
        file list from ch_bamList_f
        file '*' from ch_all_bams               
    output: 
        file("mpileup_out.vcf.gz") into ch_mpileupOUT           
    
    publishDir path: './callnocall', mode: 'copy'                                    
    
    module      samtoolsModule
    module      bcftoolsModule
    
    script:

    """
    bcftools mpileup --threads ${task.cpus} -Oz -d 250 -B -R ${exome_bed} -f ${ref} -b ${list} -o mpileup_out.vcf.gz
    """
}

//stupid un-needed GATK section
process generateBqsrModel {

    label 'gatk_unknown'

    input:
        set baseName, file(sortedBam), file(bamIndex) from ch_forGATK
    output:
        set baseName, file(sortedBam), file(bamIndex), file("${baseName}.recalreport") into ch_recalReportsBams

    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T BaseRecalibrator \
        -R $ref \
        -I $sortedBam \
        -o ${baseName}.recalreport \
        -knownSites $millsIndels \
        -knownSites $knownIndels \
        -knownSites $dbSNP \
        -L ${panel_int}
    """
}

process applyBqsrModel {

    label 'gatk_unknown'

    input:
        set baseName, file(sortedBam), file(bamIndex), file(recalReport) from ch_recalReportsBams
    output:
        set baseName, file("${baseName}.mapped.marked.sorted.recal.bam") into ch_recalibratedBams

    module      gatkModule

    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T PrintReads -R ${ref} -I ${sortedBam} --BQSR ${recalReport} \
                   -o ${baseName}.mapped.marked.sorted.recal.bam
    """
}

process call_variants{

    label 'gatk_unknown'

    input:
        set baseName, file(bam) from ch_recalibratedBams
    output:
        set baseName, file("${baseName}.mapped.marked.sorted.recal.g.vcf") into ch_gVcfs
    
    publishDir path: './output/variants/GATK/gvcf', mode: 'copy'

    module      gatkModule

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -jar $gatkJar -T HaplotypeCaller -R ${ref} --min_base_quality_score 20 \
                    --emitRefConfidence GVCF \
                    --num_cpu_threads_per_data_thread 8 \
                    -I ${bam} -L ${padded_int} -o "${baseName}.mapped.marked.sorted.recal.g.vcf"
    """
}

process mergeGVCFS {
    
    label 'gatk_unknown'

    input:
        set baseName, file(vcf) from ch_gVcfs
    output:
        set baseName, file("combined.g.vcf") into ch_combinedGVCF

    module gatkModule

    script:
    myfiles = vcf.collect().join('-V ')
    """
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T CombineGVCFs -R ${ref} \
                    --disable_auto_index_creation_and_locking_when_reading_rods \
                    -V $myfiles -o "combined.g.vcf"
    """
}

process genotypeGVCF {
    
    label 'gatk_unknown'
    
    input:
        file(vcf) from ch_combinedGVCF
    output:
        file("genotyped.vcf") into ch_genotypedGVCFsnp,ch_genotypedGVCFindel, ch_forCombining
    
    module gatkModule

    script:
    """
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T GenotypeGVCFs -R ${ref} \
                    --disable_auto_index_creation_and_locking_when_reading_rods \
                    --variant $vcf --out "genotyped.vcf" \
                    --variant $ceu_mergeGvcf --variant $gbr_mergeGvcf --variant $fin_mergeGvcf
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
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T VariantRecalibrator \
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
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods\
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
    
    script:
    """   
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T ApplyRecalibration \
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
        file("indel_recal.vcf") into ch_indeRecal    
    
    script:
    """
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods \
                    -R $ref --ts_filter_level 99.0 --excludeFiltered --num_threads ${task.cpus} \
                    -input $vcf -recalFile $recal -tranchesFile $tranches \
                    -mode INDEL -o "indel_recal.vcf"

    """
}

process combineAllRecal {
    
    label 'gatk_unknown'

    input:
        file(snp_recal) from ch_applysnpRecal
        file(indel_recal) from ch_applyindelRecal
    output:
        file("recalibrated.bam") into ch_finalGATKvcf

    script:
    """
    java -jar $gatkJar -Xmx${task.memory.toGiga() - 2}g -T CombineVariants \
                    -R $ref --disable_auto_index_creation_and_locking_when_reading_rods \
                    --num_threads ${task.cpus} --genotypemergeoption UNSORTED --variant $snp_recal \
                    --variant $indel_recal -o "recalibrated.bam"
    """
}

process chunkBEDfile {
    
    label "small_1"

    output: 
    file("*.bed") into ch_bedSegments

    module 'bedtools/2.27.1-gcc5'

    script:
    """
    bedtools split -i "$padded_bed" -n 12
    """
}

//ch_bedSegments = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")

ch_vardict= ch_forVARDICT.combine(ch_bedSegments)

process vardict {
    
    label 'medium_6h'

    input:
        set baseName, file(bam), file(bai), file(segment) from ch_vardict
    output:
        set baseName, file("${baseName}.${segment}.vardict.tsv") into ch_vardictSegments

    script:
    """
    export PATH=/home/jste0021/scripts/VarDict-1.5.8/bin/:$PATH
    VarDict -G ${ref} -f 0.01 -N "$baseName" \
        -b "$bam" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${baseName}.${segment}.vardict.tsv"
    """
}

ch_collatedSegments = ch_vardictSegments.map{ sample, segment -> [sample, segment]}.groupTuple(by: [0])

process catSegments {
    
    label 'small_short'
    echo true

    input: 
        set baseName, file(tsv) from ch_collatedSegments
    output: 
        set baseName, file("${baseName}.collated.vardict.tsv") into ch_vardictCollated
    
    script:
    myfiles = tsv.collect().join(' ')
    """
    cat ${myfiles} > ${baseName}.collated.vardict.tsv
    """
}

process makeVCF {
    
    label "big_6h"

    input:
        set sample, file(tsv) from ch_vardictCollated
    output:
        set sample, file("${sample}.vardict.vcf") into ch_outputVCF
    
    publishDir path: './output/variants/vardict', mode: 'copy'
    
    script:
    """  
    module purge
    module load R/3.5.1
    cat $tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/teststrandbias.R | \
    /home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_valid.pl -N "$sample" -f 0.01 > "${sample}.vardict.vcf"
    """
}

process collectHSMetrics {

    label 'medium_6h'

    input:
        set sample, file(bam), file(bai) from ch_forHSMetrics
    output:
        set sample, file("*.HSmetrics.txt"), file("*.perbase.txt"), file("*.pertarget.txt") into ch_metrics_unused2
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    script:
    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${sample}.HSmetrics.txt" \
        -R ${ref} \
        -BI $panel_int \
        -TI $padded_int \
        --PER_BASE_COVERAGE "${sample}.perbase.txt" \
        --PER_TARGET_COVERAGE "${sample}.pertarget.txt"
    """
}
