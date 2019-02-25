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
    for t in sequence_tuple_list[1:]:
        if temp_size + t[1] <= longest_sequence:
            temp_size += t[1]
            tsv_string += '\t' + t[0] + hg38_protection_tag
        else:
            with open("{}.contig_group".format(group_count), 'w') as out_fh:
                out_fh.write(tsv_string)
            group_count += 1
            tsv_string = t[0] + hg38_protection_tag
            temp_size = t[1]
    if tsv_string:
        with open("{}.contig_group".format(group_count), 'w') as out_fh:
            out_fh.write(tsv_string)
        group_count += 1
    with open("{}.contig_group".format(group_count), 'w') as out_fh:
        out_fh.write("unmapped")
    """
}


// See https://github.com/nextflow-io/nextflow/issues/298 for combine 
// operator to get cartesian product of two channels
contigBamScatter_ch = sortedBamFiles.combine(contigGroupings.flatten())


/* process printInputs { */
/*     stageInMode 'symlink' */

/*     input: */
/*         set baseName, sortedBam, bamIndex, contigGrouping  from contigBamScatter_ch */
/*     output: */
/*         stdout result */

/*     """ */
/*     # Dirty way of reading and formatting the contig list */
/*     CONTIG_LIST=\$(cat $contigGrouping | sed -E 's/\t|^/ INPUT=/g') */
/*     echo \$CONTIG_LIST */
/*     """ */
/* } */


process generateBqsrModel {
}

// These could possibly be merged? Would mean for a longer job though

process applyBqsrModel {
}


/* process gatherBqsrReports { */
/* } */


/* process gatherBams { */
/* } */


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

// "${refFolder}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
// "${refFolder}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
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
