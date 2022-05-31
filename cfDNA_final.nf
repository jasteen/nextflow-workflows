#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./test')
panel_int      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions_b37.interval_list')
padded_int     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded_b37.interval_list')
panel_bed      = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Regions_b37.bed')
padded_bed     = file('/projects/vh83/reference/sureselect/medha_exome_panel/S30409818_Padded_b37.bed')
tmp_dir        = file('/scratch/vh83/tmp/')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

//VEP
//Annotation resources
vep_cache      = file("/projects/vh83/reference/VEP_CACHE")


// Tools
picardJar          = '~/picard.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
gatkModule         = 'gatk/4.0.11.0' 
rModule            = 'R/3.5.1'          
fgbioJar           = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'
condaModule        = 'miniconda3/4.1.11-python3.5' 

// Creating channel from input directory
Channel.fromFilePairs("$inputDirectory/*_{R1,R2}.fastq.gz", size: 2, flat: true).into{ch_inputFiles;ch_forFastqc}

process runFASTQC {
    label 'small_2'

    input:
        set baseName, file(R1), file(R2) from ch_forFastqc
    output:
        file("*.{html,zip}") into ch_fastqcReports
    
    publishDir path: './output/metrics/fastqc', mode: 'copy'

    module      'fastqc'
        
    script:
    
    """
    fastqc -t ${task.cpus} -q $R1 $R2
    """

}

process surecallTrimmer {
    
    label 'genomics_1'

    input:
        set baseName, file(R1) file(R2) from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped.umi.bam") into ch_unmappedUMIbams

    module 'java/openjdk-1.14.02' 
    """
    bash /projects/vh83/local_software/agent3.0/agent.sh trim -fq1 ${R1} -fq2 ${R2} -v2 -bam -out ./${baseName}.unmapped.umi  
    """
       
}


process alignBwa {
    
    label 'bwa'
    
    input:
        set baseName, file(bam) from ch_unmappedUMIbams
    output:
        set baseName, file("${baseName}.aligned.bam") into ch_pipedBams, ch_forMetrics1

    publishDir path: './output/bams', mode: 'copy'

    module      bwaModule
    module	    'samtools'
    module      'picard'
    
    script:
    """
    set -o pipefail
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ '/dev/stdout' -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 \
        -INTERLEAVE true -NON_PF true -TMP_DIR "$tmp_dir" | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx${(task.memory.toGiga() / 6).toInteger() }g -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.aligned.bam" -R "$ref" -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process groupreadsByUmi {
    
    label 'medium_6h'

    input:
        set baseName, file(bam) from ch_pipedBams
    output:
        set baseName, file("${baseName}.piped.grouped.histogram.tsv"), file("${baseName}.piped.grouped.bam") into ch_umiGroupedBams
    
    publishDir path: './output/metrics/UMI/family_sizes', mode: 'copy', pattern: "*.tsv"

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar GroupReadsByUmi \
         -i ${bam} -f "${baseName}.piped.grouped.histogram.tsv" -o "${baseName}.piped.grouped.bam" -s Adjacency -e 1 
    """

}

process generateConsensusReads {
    
    label 'medium_6h'

    input:
        set baseName, file(hist), file(bam) from ch_umiGroupedBams
    output:
        set baseName, file("${baseName}.consensus.unmapped.bam") into ch_unmappedConsensusBams
    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar CallMolecularConsensusReads \
        --input $bam --output ${baseName}.consensus.unmapped.bam \
        --error-rate-post-umi 30 --min-reads 3
    """
}

process mapConsensusReads {
    
    label 'bwa'

    input:
        set baseName, file(bam) from ch_unmappedConsensusBams
    output:
        set baseName, file("${baseName}.consensus.aligned.bam") into ch_mappedConsensusBams, ch_forMetrics2
    publishDir path: './output/bams', mode: 'copy'

    module 	    bwaModule
    
    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR $tmp_dir | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.consensus.aligned.bam" -R $ref -ADD_MATE_CIGAR true \
        -SO coordinate -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process indexBam {
    
    label 'small_1'

    input:
        set baseName, file(bam) from ch_mappedConsensusBams
    output:
        set baseName, file(bam), file("${baseName}.consensus.aligned.bam.bai") into ch_indexedConsensusBams
    publishDir path: './output/bams', mode: 'copy'
 
    module      'samtools'
 
    script:
    """
    samtools index $bam ${baseName}.consensus.aligned.bam.bai
    """

}

//I cant think of a better workflow to make sure Tumor and Normal are processed in order
//create channel for normal and tumor
ch_tumor  = Channel.create()
ch_normal = Channel.create()

//split single bam channel into tumor and normal **CURRENTLY RELIES ON "SAMPLE_[FFPE|NORMAL]" naming scheme
ch_indexedConsensusBams.choice(ch_tumor, ch_normal){ a -> a[0] =~ /cfDNA/ ? 0 : 1 }

//split SAMPLE from FFPE|NORMAL so channels can be joined by sample
ch_normalSplit = ch_normal.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}
ch_tumorSplit = ch_tumor.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}

//merge tumor and normal back together by sample number. 
ch_tumorNormalPairs = ch_tumorSplit.join(ch_normalSplit)

//create bedfile segments
ch_bedSegments = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")
//create cartesian product of the input channel and the segments files
ch_vardictInput = ch_tumorNormalPairs.combine(ch_bedSegments)

process runVardict {
    
    label 'vardict'

    input:
        set sample, ttype, file(tbam), file(tbai), ntype, file(nbam), file(nbai), file(segment) from ch_vardictInput
    output:
        set sample, file(tbam), file(nbam), file("${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv") into ch_rawVardictSegments
 
    script:
    """
    export PATH=/home/jste0021/scripts/VarDict-1.7.0/bin/:$PATH
    VarDict -G ${ref} -f 0.01 -N "${tbam}|${nbam}" \
        -b "${tbam}|${nbam}" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv"
    """ 

}

ch_collatedSegments = ch_rawVardictSegments.map{ sample, tbam, nbam, segment -> [sample, tbam.name, nbam.name, segment]}.groupTuple(by: [0,1,2])

process catSegments {
    
    label 'small_1'
    echo true

    input: 
        set sample, tbam, nbam, file(tsv) from ch_collatedSegments
    output: 
        set sample, tbam, nbam, file("${sample}.collated.vardict.tsv") into ch_rawVardict

    script:  
    myfiles = tsv.collect().join(' ')
    """
    cat ${myfiles} > ${sample}.collated.vardict.tsv
    """

}

process makeVCF {
    
    label 'medium_6h'

    input:
        set sample, tbam, nbam, file(tsv) from ch_rawVardict
    output:
        set sample, file("${sample}.somatic.vardict.vcf") into ch_outputVCF
    
    //publishDir path: './output/UMI/intermediate', mode: 'copy'
    
    script:
 
    """  
    module purge
    module load R/3.5.1
    cat $tsv | /home/jste0021/scripts/VarDict-1.7.0/bin/testsomatic.R | \
    /home/jste0021/scripts/VarDict-1.7.0/bin/var2vcf_paired.pl -N "${tbam}|${nbam}" -f 0.00001 > "${sample}.somatic.vardict.vcf"
    """
}

process reheaderUMIVCF {
    
    label 'small_1'

    input:
        set sample, file(vcf) from ch_outputVCF
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

process sortVCFS {

    label 'medium_6h'
    
    input:
        set baseName, file(vcf) from ch_reheaderVCF
    output:
        set baseName, file("${baseName}.UMI.reheader.sorted.vcf.gz") into ch_sortedVCF

    //publishDir path: './output/vcf/UMI', mode: 'copy'                                    
    
    module     'bcftools/1.8'                       
    
    script:
    """
    bcftools sort -o "${baseName}.UMI.reheader.sorted.vcf.gz" -O z ${vcf}
    """
}

process indexVCFS {

    label 'small_1'

    input:
        set baseName, file(vcf) from ch_sortedVCF
    output:
        set baseName, file(vcf), file("${baseName}.UMI.reheader.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './output/vcf/UMI', mode: 'copy'                                    
    
    module     'bcftools/1.8'

    script:
  
    """
    bcftools index -f --tbi ${vcf} -o ${baseName}.UMI.reheader.sorted.vcf.gz.tbi
    """
}

process vt_decompose_normalise {
        
    label 'medium_6h'

    input:
        set baseName, file(vcf), file(tbi) from ch_indexedVCF
    output:
        set baseName, file("${baseName}.UMI.reheader.sorted.vt.vcf.gz") into ch_vtDecomposeVCF

    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    module      'vt'
    
    script:
    """
    vt decompose -s $vcf | vt normalize -r $ref -o "${baseName}.UMI.reheader.sorted.vt.vcf.gz" -
    """
}

process apply_vep {
    label 'vep_sing'

    input:
        set baseName, file(vcf) from ch_vtDecomposeVCF
    output:
        set baseName, file("${baseName}.UMI.VEP_Stats.html"), file("${baseName}.UMI.reheader.sorted.vt.vep.vcf") into ch_vepVCF

    publishDir path: './output/vcf/UMI', mode: 'copy', pattern: "*.vcf"
    publishDir path: './output/metrics/vep_stats', mode: 'copy', pattern: "*.html"

    module 'singularity'

    script:
    """
    vep --cache --dir_cache $vep_cache \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --stats_file "${baseName}.UMI.VEP_Stats.html"  \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o "${baseName}.UMI.reheader.sorted.vt.vep.vcf"
    """
}



//ch_forMetrics = ch_forMetrics1.concat(ch_forMetrics2)
ch_forMetrics1.concat(ch_forMetrics2).into{ch_forMultipleMetrics;ch_forHSMetrics}

process collectHSMetrics {

    label 'medium_6h'

    input:
        set sample, file(bam) from ch_forHSMetrics
    output:
        set sample, file("*.HSmetrics.txt"), file("*.pertarget.txt") into ch_metrics
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${bam.baseName}.HSmetrics.txt" \
        -R ${ref} \
        -BI $panel_int \
        -TI $padded_int \
        --PER_TARGET_COVERAGE "${bam.baseName}.pertarget.txt"
    """
}

process collectMultipleMetrics {

    label 'medium_6h'

    input:
        set sample, file(bam) from ch_forMultipleMetrics
    output:
        set sample, file("*multiple_metrics*") into ch_metrics2
    
    publishDir path: './output/metrics/multiple', mode: 'copy'
    
    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx${task.memory.toGiga() - 2}g -jar ${picardJar} CollectMultipleMetrics \
        -I $bam \
        -O ${bam.baseName}.multiple_metrics \
        -R $ref
    """
}

process multiQC {

    label 'medium_6h'

    input:
        file('coverage/*') from ch_metrics2.collect()
        file('multiple/*') from ch_metrics.collect()
        file('fastqc/*') from ch_fastqcReports.collect()
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
