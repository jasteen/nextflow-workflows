#!/usr/bin/env nextflow

// Reference Files
refFolder = file("/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/")
refBase   = "$refFolder/Homo_sapiens_assembly38"
ref       = file("${refBase}.fasta")
refDict   = file("${refBase}.dict")
refFai    = file("${refBase}.fasta.fai")

// Tools
// Note: Necessary to provide absolute path to jar files to set memory limits
//       in some cases.
picardJar = "/usr/local/picard/2.9.2/bin/picard.jar"

// Inputs
inputDirectory = file('/scratch/vh83/sandbox/jared/full_cwl_pipeline_testing/input_files/')
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz").take(1)

// Outputs
outputDir = "/scratch/vh83/sandbox/jared/full_cwl_pipeline_testing/nextflow/outputs"


process alignBwa {
    stageInMode 'symlink'

    input:
        set baseName, file(fastqs) from inputFiles
    output:
        set baseName, file("${baseName}.bam") into bamFiles

    module 'bwa/0.7.17-gcc5'
    module 'samtools/1.9'

    // TODO: This should result in queryname sorted output but isn't for some
    //       reason. Could be a version issue with samtools or bwa.
    //       Replace sort with "samtools view -b -h -o ${baseName}.bam -" if
    //       fixed.
    """
    set -o pipefail
    bwa mem \
        -K 100000000 -v 3 -Y -t 1 \
        -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" \
        $ref ${fastqs[0]} ${fastqs[1]} | \
        java -Xms4000m -jar $picardJar SortSam \
            INPUT=/dev/stdin \
            OUTPUT=${baseName}.bam \
            SORT_ORDER=queryname
    """
}


process markDuplicatesPicard {
    stageInMode 'symlink'
    publishDir outputDir

    input:
        set baseName, bam from bamFiles 
    output:
        set baseName, file("${baseName}.marked.bam") into markedBamFiles
        set baseName, file("${baseName}.markduplicates.metrics") into metrics

    // TODO: CLEAR_DT=false option in GATK pipeline but not supported by 
    //       this version of picard.
    //       ADD_PG_TAG_TO_READS=false also not supported.
    module  'picard/2.9.2'
    """
    java -Xms4000m -jar $picardJar MarkDuplicates \
        INPUT=$bam \
        OUTPUT=${baseName}.marked.bam \
        METRICS_FILE=${baseName}.markduplicates.metrics \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER=queryname
    """
}




// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.alt"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.amb"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.sa"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.pac"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.dict"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.bwt"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.fai"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.fasta.ann"

// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/hapmap_3.3.hg38.vcf.gz"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/wgs_calling_regions.hg38.interval_list"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G_omni2.5.hg38.vcf.gz.tbi" 
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/1000G_omni2.5.hg38.vcf.gz"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/hapmap_3.3.hg38.vcf.gz.tbi"
// "/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
