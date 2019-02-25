#!/usr/bin/env nextflow

// Reference Files
refFolder = file("/projects/vh83/reference/genomes/hg38/hg38_broad_resource_bundle/v0/")
refBase   = "$refFolder/Homo_sapiens_assembly38"
ref       = file("${refBase}.fasta")
refDict   = file("${refBase}.dict")
refFai    = file("${refBase}.fasta.fai")

// Known Sites
millsIndels = file("${refFolder}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
knownIndels = file("${refFolder}/Homo_sapiens_assembly38.known_indels.vcf.gz")
dbSNP       = file("${refFolder}/Homo_sapiens_assembly38.dbsnp138.vcf")


// Tools
// Note: It is necessary to provide absolute path to jar files to set memory limits
//       in some cases.
// TODO: Set file versions here to be used by processes (allows us to update 
//       versions easily).
picardJar = "/usr/local/picard/2.9.2/bin/picard.jar"


// Inputs
inputDirectory = file('/scratch/vh83/sandbox/jared/full_cwl_pipeline_testing/input_files/')
inputFiles = Channel.fromFilePairs("$inputDirectory/*_R{1,2}.fastq.gz").take(2)

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


process sortBam {
    stageInMode 'symlink'

    input:
        set baseName, markedBam from markedBamFiles
    output:
        set baseName,
            file("${baseName}.marked.sorted.bam"), 
            file("${baseName}.marked.sorted.bai") into sortedBamFiles

    module  'picard/2.9.2'
    """
    java -Xms4000m -jar $picardJar SortSam \
        INPUT=$markedBam \
        OUTPUT=${baseName}.marked.sorted.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000
    """
}


process groupContigs {
    output:
        file("*.contig_group") into contigGroupings

    // Reads the dict file and divides into tab separated contig groupings
    // One per file
    // Adapted from CreateSequenceGroupingTSV at: 
    // https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.gatk4.0.wdl
    // TODO: Clean up & rewrite where appropriate
    module 'python/3.7.2-gcc6'
    """
    #!/usr/bin/env python
    with open("$refDict", 'r') as fh:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in fh:
            if line.startswith("@SQ"):
                sl = line.split('\t')
                sequence_tuple_list.append((sl[1].split("SN:")[1],
                                            int(sl[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # GATK4 strips off anything following the final ':', adding this to preserve
    # contigs containing ':'
    hg38_protection_tag = ":1+"
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    group_count = 0

    def write_group(tsv_string, group_count):
        group_count_str = str(group_count)
        n_zeros = 3 - len(group_count_str)
        group_count_str = (n_zeros * "0") + group_count_str
        with open("{}.contig_group".format(group_count_str), 'w') as out_fh:
            out_fh.write(tsv_string)


    for t in sequence_tuple_list[1:]:
        if temp_size + t[1] <= longest_sequence:
            temp_size += t[1]
            tsv_string += '\t' + t[0] + hg38_protection_tag
        else:
            write_group(tsv_string, group_count)
            group_count += 1
            tsv_string = t[0] + hg38_protection_tag
            temp_size = t[1]
    if tsv_string:
        write_group(tsv_string, group_count)
        group_count += 1
    write_group("unmapped", group_count)
    """
}


// Using combine operator to get cartesian product of two channels
contigBamScatter_ch = sortedBamFiles.combine(contigGroupings.flatten())


process generateBqsrModel {
    input:
        set baseName, sortedBam, bamIndex, contigGrouping  from contigBamScatter_ch
    output:
        set baseName, sortedBam, bamIndex, contigGrouping, 
            file("${baseName}.${contigGrouping.baseName}.recalreport") into recalReportsBams

    module 'gatk/4.0.11.0'
    """
    # Somewhat dirty way of formatting the sequence groupings
    SEQUENCE_GROUPING=\$(cat $contigGrouping | sed -E 's/\t|^/ -L /g')
    gatk --java-options '-Xms4000m' BaseRecalibrator \
        --use-original-qualities \
        -R $ref \
        -I $sortedBam \
        -O ${baseName}.${contigGrouping.baseName}.recalreport \
        --known-sites $millsIndels \
        --known-sites $knownIndels \
        --known-sites $dbSNP \
        \$SEQUENCE_GROUPING
    """
}


process applyBqsrModel {
    input:
        set baseName, sortedBam, bamIndex, contigGrouping, recalReport from recalReportsBams
    output:
        set baseName, contigGrouping,
            file("${baseName}.${contigGrouping.baseName}.recal.bam") into recalibratedBams

    module 'gatk/4.0.11.0'
    """
    # Somewhat dirty way of formatting the sequence groupings
    SEQUENCE_GROUPING=\$(cat $contigGrouping | sed -E 's/\t|^/ -L /g')
    gatk --java-options '-Xms3000m' ApplyBQSR \
    --add-output-sam-program-record \
    --use-original-qualities \
    --static-quantized-quals 10 \
    --static-quantized-quals 20 \
    --static-quantized-quals 30 \
    -R $ref \
    -I $sortedBam \
    -O ${baseName}.${contigGrouping.baseName}.recal.bam \
    -bqsr $recalReport \
    \$SEQUENCE_GROUPING
    """
}


// Create a new channel grouping the recalibrated fragments by baseName (first
// part of tuple)
groupedRecalibratedBams_ch = recalibratedBams.groupTuple()



process gatherBams {
    input:
        set baseName, contigGroupings, recalBams from groupedRecalibratedBams_ch
    output:
        set baseName, file("${baseName}.recal.merge.bam") into recalMergedBams

    module  'picard/2.9.2'
    // Note: Due to the use of temp folders the recalibrated fragments do not
    //       naturally sort correctly. As order is important to retain sorting
    //       they are sorted by their basename prior to being added as inputs.
    // TODO: Do something with metrics file. Or update picard. Seems newer
    //       versions don't require it as an input.
    """
    java -Xms4000m -jar $picardJar MarkDuplicates \
        OUTPUT=${baseName}.recal.merge.bam \
        CREATE_INDEX=true \
        METRICS_FILE=metrics.txt \
        ${" INPUT=" + recalBams.sort { it.baseName }.join(" INPUT=")}
    """
}



/* process calculateScatterIntervals { */
/* } */
/* // Cartesian product of these */
/* process callHaplotypeCallerGvcf { */
/* } */


/* // Merge haplotypecaller output keyed by file baseName */
/* process mergeGvcfs { */
/* } */



// Many more metrics processes to be added but this is fine for now





// "${refFolder}/Homo_sapiens_assembly38.fasta"
// "${refFolder}/Homo_sapiens_assembly38.fasta.alt"
// "${refFolder}/Homo_sapiens_assembly38.fasta.amb"
// "${refFolder}/Homo_sapiens_assembly38.fasta.sa"
// "${refFolder}/Homo_sapiens_assembly38.fasta.pac"
// "${refFolder}/Homo_sapiens_assembly38.dict"
// "${refFolder}/Homo_sapiens_assembly38.fasta.bwt"
// "${refFolder}/Homo_sapiens_assembly38.fasta.fai"
// "${refFolder}/Homo_sapiens_assembly38.fasta.ann"

// "${refFolder}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
// "${refFolder}/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx"
// "${refFolder}/Homo_sapiens_assembly38.dbsnp138.vcf"
// "${refFolder}/hapmap_3.3.hg38.vcf.gz"
// "${refFolder}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
// "${refFolder}/wgs_calling_regions.hg38.interval_list"
// "${refFolder}/1000G_omni2.5.hg38.vcf.gz.tbi" 
// "${refFolder}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
// "${refFolder}/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"
// "${refFolder}/Homo_sapiens_assembly38.known_indels.vcf.gz"
// "${refFolder}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
// "${refFolder}/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
// "${refFolder}/1000G_omni2.5.hg38.vcf.gz"
// "${refFolder}/hapmap_3.3.hg38.vcf.gz.tbi"
// "${refFolder}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
