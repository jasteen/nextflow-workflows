#!/usr/bin/env nextflow

// Required Inputs
inputDirectory = file('./fastqs/')
params.reference = ""

params.panel_bed       = ""
params.padded_bed      = ""
params.panel_int       = ""
params.padded_int      = ""

//set up for multiple regerence possibilities.

if(params.reference == "hg19"){
    //HG19 reference for aspree stuff
    refFolder      = file("/projects/vh83/reference/genomes/hg19/bwa_0.7.12_index")
    refBase        = "$refFolder/ucsc.hg19"
    ref            = file("${refBase}.fasta")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")

}else if(params.reference == "hg38"){
    refFolder      = file("/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0")
    refBase        = "$refFolder/Homo_sapiens_assembly38"
    ref            = file("${refBase}.fasta")
    refDict        = file("${refBase}.dict")
    refFai         = file("${refBase}.fasta.fai")
    genome_file    = file("$refFolder/hg38.chrom.sizes")
    header         = file("$refFolder/hg38_vcf_header.txt")
    vep_cache      = file("/projects/vh83/reference/VEP_CACHE")


}else{
    
    refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index")
    refBase          = "$refFolder/human_g1k_v37_decoy"
    ref              = file("${refBase}.fasta")
    refDict          = file("${refBase}.dict")
    refFai           = file("${refBase}.fasta.fai")
}



// Tools
picardJar      = '/usr/local/picard/2.19.0/bin/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'
surecalltrimmerJar = '/projects/vh83/local_software/agent/SurecallTrimmer_v4.0.1.jar'
condaModule        = 'miniconda3/4.1.11-python3.5'
rModule            = 'R/3.5.1'
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

// Creating channel from input directory
Channel.fromFilePairs("$inputDirectory/*_R{1,2}_001.fastq.gz") into {ch_processedInputFiles;ch_forFastqc}


process fastQC {
    
    label 'start_1_8_2h'

    publishDir path: './fastqc_results', mode: 'copy'
    
    input:
        set baseName, file(fastq) from ch_forFastqc
    output:
        set baseName, file("${baseName}*.zip"), file("${baseName}*.html") into ch_fastQCOutput

    module      'fastqc'
    
    """
    fastqc $fastq
    """
}


process align_bwa {

    label 'bwa'

    input:
        set baseName, file(fastqs) from ch_processedInputFiles
    output:
        set baseName, file("${baseName}.hq.sorted.bam"), file("${baseName}.hq.sorted.bam.bai") into ch_mappedBams

    publishDir path: './bam_out', mode: 'copy'

    module      bwaModule
    module      samtoolsModule

    """
    bwa mem -M -t ${task.cpus} -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" $ref ${fastqs[0]} ${fastqs[1]}  \
        | samtools view -u -h -q 1 -f 2 -F 4 -F 8 -F 256 - \
        | samtools sort -@ ${task.cpus} -o "${baseName}.hq.sorted.bam"
    samtools index "${baseName}.hq.sorted.bam" "${baseName}.hq.sorted.bam.bai"
    """
}

process markDuplicatesPicard {
    
    label 'start_1_8_15m'

    input:
        set baseName, bam, bai from ch_mappedBams
    output:
        set baseName, file("${baseName}.hq.sorted.marked.bam"), file("${baseName}.hq.sorted.marked.bam.bai") into ch_markedBamFiles
        set baseName, file("${baseName}.markduplicates.metrics") into ch_metrics_unused

    publishDir path: './output/metrics/markduplicates', mode: 'copy', pattern: '*.metrics'
    publishDir path: './output/bams/markdup', mode: 'copy', pattern: '*.{bam,bai}'


    module      samtoolsModule
    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar $picardJar MarkDuplicates \
        -INPUT $bam \
        -OUTPUT ${baseName}.hq.sorted.marked.bam \
        -METRICS_FILE ${baseName}.markduplicates.metrics \
        -VALIDATION_STRINGENCY SILENT \
        -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        -ASSUME_SORT_ORDER queryname
    samtools index "${baseName}.hq.sorted.marked.bam" "${baseName}.hq.sorted.marked.bam.bai"

    """
}

ch_markedBamFiles.into{ch_mappedBam1;ch_mappedBam2;ch_mappedBam3;ch_mappedBam4;ch_mappedBam5;ch_mappedBam6;ch_forMetrics1}


process run_vardict {

    label 'start_1_8_2h'

    input:
        set baseName, file(bam), file(bai) from ch_mappedBam1               
    output: 
        set baseName, file("${baseName}.tsv") into ch_vardict_TSV           
    
    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    script:
    """
    export PATH=/home/jste0021/scripts/git_controlled/vardict_testing/VarDictJava/build/install/VarDict/bin/:$PATH
    VarDict -G ${ref} -f 0.1 -N "${baseName}" --nosv -b ${bam} -c 1 -S 2 -E 3 -g 4 ${params.panel_bed} > "${baseName}.tsv"
    """
}

process makeVCF {

    label 'start_1_8_2h'

    input:
        set baseName, file(tsv) from ch_vardict_TSV
    output:
        set baseName, file("${baseName}.vardict.vcf") into ch_vardictVCFs
    
    publishDir path: './variants_raw_out', mode: 'copy'
    
    script:
    """
    module purge
    module load R/3.5.1
    cat ${tsv} | /home/jste0021/scripts/VarDict-1.7.0/bin/teststrandbias.R | \
        /home/jste0021/scripts/VarDict-1.7.0/bin/var2vcf_valid.pl -N "${baseName}" \
        -f 0.1 -E > "${baseName}.vardict.vcf"
    """
}

process reheaderVCF {

    label 'start_1_8_2h'

    input:
        set baseName, file(vcf) from ch_vardictVCFs
    
    output:
        set baseName, file("${baseName}.reheader.vcf.gz") into ch_reheaderVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
    """
    bcftools annotate -h ${header} -O z -o "${baseName}.reheader.vcf.gz" ${vcf}
    """

}

process sortVCFS {

    label 'start_1_8_2h'

    input:
        set baseName, file(vcf) from ch_reheaderVCF
    output:
        set baseName, file("${baseName}.sorted.vcf.gz") into ch_sortedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
    """
    bcftools view ${vcf} | awk 'BEGIN{FS=OFS="\t"}{if(\$0 ~ /^#/) print;else if(\$4 ~ /[ACTG]/) print}' | bcftools view -O z -o ${baseName}.fixed.vcf.gz
    bcftools sort -o "${baseName}.sorted.vcf.gz" -O z ${baseName}.fixed.vcf.gz
    """
}

process indexVCFS {

    label 'start_1_8_2h'

    input:
        set baseName, file(vcf) from ch_sortedVCF
    output:
        set baseName, file(vcf), file("${baseName}.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
    """
    bcftools index -f --tbi ${vcf} -o ${baseName}.sorted.vcf.gz.tbi
    """
}

//duplicate ch_indexedVCF
ch_indexedVCF.into{ch_list;ch_files}
//set one version to a list of filenames of the VCF
ch_list.map {it -> it[1].name}
    .collectFile(name: 'list2.txt', newLine: true)
    .set {ch_list_f}
//set the second to all the files
ch_files
    .collect()
    .set {ch_all_files}

//feed both to the merge so that the indexes are available to bcftools

process mergeVCFS {

    label 'start_1_8_2h'

    echo true

    publishDir './variants_merged/', mode: 'copy'

    input:
    file list from ch_list_f
    file '*' from ch_all_files
    
    output:
    file "merged_vardict_norm.vcf.gz" into ch_mergedfinalVCF

    module     'bcftools/1.8'
    
    script: 
    """
    split -l 500 list2.txt temp_shorter_list_
    for i in temp_shorter_*; do bcftools merge -m all -l \$i -O z -o \$i.merged.vcf.gz; bcftools index -t \$i.merged.vcf.gz; done
    ls *merged.vcf.gz > list3.txt
    bcftools merge -m all -O z -o "merged_vardict.vcf.gz" -l list3.txt
    bcftools norm -f $ref -m -both -o "merged_vardict_norm.vcf.gz" -O z "merged_vardict.vcf.gz"
    bcftools index --tbi "merged_vardict_norm.vcf.gz"
    """
}


process apply_vep {
    label 'vep_sing'

    input:
        file(vcf) from ch_mergedfinalVCF
    output:
        file("merged_vardict_norm_vep.vcf") into ch_vepVCF

    publishDir path: './variants_merged', mode: 'copy'

    module 'singularity'

    script:
    """
    vep --cache --dir_cache $vep_cache \
                      --assembly GRCh38 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o merged_vardict_norm_vep.vcf
    """
}
    

/*
Stats Generation Section
*/

ch_forMetrics1.into{ch_forMultipleMetrics;ch_forHSMetrics}

process collectHSMetrics {

    label 'start_1_8_2h'

    input:
        set sample, file(bam), file(bai) from ch_forHSMetrics
    output:
        set sample, file("*.HSmetrics.txt"), file("*.pertarget.txt") into ch_metrics
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    script:

    """
    module purge
    module load R/3.5.1
    module picard/2.19.0

    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${bam.baseName}.HSmetrics.txt" \
        -R ${ref} \
        -BI ${params.panel_int} \
        -TI ${params.padded_int} \
        --PER_TARGET_COVERAGE "${bam.baseName}.pertarget.txt"
    """
}

process multiQC {

    label 'start_1_8_2h'

    input:
        file('coverage/*') from ch_metrics.collect()
        file('fastqc/*') from ch_fastQCOutput.collect()
        
    output:
        set file("*multiqc_report.html"), file("*multiqc_data") into ch_multiQCOut

    publishDir path: './output/metrics/report', mode: 'copy'

    module      condaModule
    conda       '/home/jste0021/.conda/envs/py3.5/'

    script:
    
    """
    multiqc -f -v .
    """
}
