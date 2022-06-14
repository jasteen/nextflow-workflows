#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./test')


params.panel_bed       = ""
params.padded_bed      = ""
params.panel_int       = ""
params.padded_int      = ""

tmp_dir        = file('/scratch/vh83/tmp/')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.00001
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
fgbioJar           = '/fs02/vh83/local_software/fgbio/fgbio-2.0.2.jar'
condaModule        = 'miniconda3/4.1.11-python3.5' 

// Creating channel from input directory
Channel.fromFilePairs("$inputDirectory/*_{R1,R2}.fastq.gz", size: 2, flat: true).into{ch_inputFiles;ch_forFastqc}

process runFASTQC {
    label 'start_1_8_2h'

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
    
    label 'start_1_16_6h'

    input:
        set baseName, file(R1), file(R2) from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped_R1.fastq.gz"), file("${baseName}.unmapped_R2.fastq.gz") into ch_surecall

    module 'java/openjdk-1.14.02' 
    module 'samtools'
    """
    bash /projects/vh83/local_software/agent3.0/agent.sh trim -fq1 ${R1} -fq2 ${R2} -v2 -out ./${baseName}.unmapped
    
    """
    
}

process alignBwa {
    
    label 'bwa'
    
    input:
        set baseName, file(R1), file(R2) from ch_surecall
    output:
        set baseName, file("${baseName}.aligned.bam") into ch_pipedBams, ch_forMetrics1

    publishDir path: './output/bams', mode: 'copy'

    module      bwaModule
    module      samtoolsModule
    
    
    script:
    """
    bwa mem -C -t ${task.cpus} -R "@RG\\tID:A\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" $ref $R1 $R2 | samtools view -b - > ${baseName}.aligned.bam

    """
}


process setMateInfo {
    label 'start_1_16_12h'

    input:
        set baseName, file(bam) from ch_pipedBams
    output:
        set baseName, file("${baseName}.aligned.matefixed.bam") into ch_mateFixed

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar SortBam \
           -i $bam -o ${baseName}.aligned.sorted.bam -s Queryname
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar SetMateInformation \
          -i ${baseName}.aligned.sorted.bam -r $ref -o ${baseName}.aligned.matefixed.bam
    """
}


process groupreadsByUmi {
    
    label 'start_1_16_12h'

    input:
        set baseName, file(bam) from ch_mateFixed
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
    
    label 'start_1_16_12h'

    input:
        set baseName, file(hist), file(bam) from ch_umiGroupedBams
    output:
        set baseName, file("${baseName}.consensus.unmapped.bam") into ch_unmappedConsensusBams, ch_duplexStats
    //publishDir path: './output/UMI/intermediate', mode: 'copy'

    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar CallMolecularConsensusReads \
        --input $bam --output ${baseName}.consensus.unmapped.bam \
        --error-rate-post-umi 30 --min-reads 3
    """
}

process generateUMIstats {
    
    label 'start_1_8_2h'

    input:
        set baseName, file(bam) from ch_duplexStats
    output:
        set baseName, file("*") into ch_stats_out
    
    publishDir path: './output/metrics/UMI', mode: 'copy'

    module       rModule
    script:
    """
    java -Xmx${task.memory.toGiga() - 2}g -Djava.io.tmpdir=$tmp_dir -jar $fgbioJar CollectDuplexSeqMetrics \
        -i $bam -o ${baseName} 
        
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
    module      samtoolsModule

    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx${ (task.memory.toGiga() / 6).toInteger() }g -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ /dev/stdout \
        -INTERLEAVE true -TMP_DIR $tmp_dir | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin > ${baseName}.temp.bam
    samtools sort -o ${baseName}.consensus.aligned.bam ${baseName}.temp.bam
    """
}

process indexBam {
    
    label 'start_1_8_15m'

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
ch_bedSegments = Channel.fromPath("${params.padded_bed}").splitText( by: 50000, file: "seg")
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
    VarDict -G ${ref} -f $af_thr -N "${tbam}|${nbam}" \
        -b "${tbam}|${nbam}" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv"
    """ 

}

ch_collatedSegments = ch_rawVardictSegments.map{ sample, tbam, nbam, segment -> [sample, tbam.name, nbam.name, segment]}.groupTuple(by: [0,1,2])

process catSegments {
    
    label 'start_1_8_15m'
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
    
    label 'start_1_8_2h'

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
    
    label 'start_1_8_15m'

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

    label 'start_1_16_6h'
    
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

    label 'start_1_8_15m'

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
        
    label 'start_1_8_2h'

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
//ch_forMetrics1.concat(ch_forMetrics2).into{ch_forMultipleMetrics;ch_forHSMetrics}
ch_forMetrics1.into{ch_forMultipleMetrics;ch_forHSMetrics}

process collectHSMetrics {

    label 'start_1_8_2h'

    input:
        set sample, file(bam) from ch_forHSMetrics
    output:
        set sample, file("*.HSmetrics.txt"), file("*.pertarget.txt") into ch_metrics2
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    script:

    """
    module purge
    module load R/3.5.1
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
        file('coverage/*') from ch_metrics2.collect()
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

