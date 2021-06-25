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
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

//VEP
//Annotation resources
dbsnp_b37       = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")
other_vep       = file("/usr/local/vep/90/ensembl-vep/cache")


// Tools
picardJar          = '~/picard.jar'
bwaModule          = 'bwa/0.7.17-gcc5'
samtoolsModule     = 'samtools/1.9'
rModule            = 'R/3.5.1'          
fgbioJar           = '/usr/local/fgbio/0.9.0/target/fgbio-0.9.0-17cb5fb-SNAPSHOT.jar'

// Global Resource Configuration Options
globalCores       = 1
bwaCores	      = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'


// Creating channel from input directory
//create channel flat because we want to join it later, and the tuple makes that more annoying than I want it to be
ch_inputFiles = Channel.fromFilePairs("$inputDirectory/*_{R1,R2}.fq.gz", size: 2, flat: true)


process createUnmappedBam {
    
    label 'genomics_2'

    publishDir path: './output/intermediate', mode: 'copy'
    
    input:
        set baseName, file(R1), file(R2) from ch_inputFiles
    output:
        set baseName, file("${baseName}.unmapped.bam") into ch_unmappedbams

    //publishDir path: './output/intermediate', mode: 'copy'
    
    
    module      'fgbio'
    module      'java'
    
    
    script:
    """
    java -Xmx30g -Djava.io.tmpdir=$tmp_dir -XX:+AggressiveOpts -XX:+AggressiveHeap \
        -jar $fgbioJar FastqToBam --input $R1 $R2 --output "${baseName}.unmapped.bam" --read-structures +T +T \
        --sample "${baseName}" --read-group-id "${baseName}" --library A --platform illumina --sort true
    """
}

process markAdaptors {

    label 'genomics_3'
    
    publishDir path: './output/metrics/UMI/adaptor_marking', mode: 'copy', pattern: '*.tsv'

    input:
        set baseName, file(bam) from ch_unmappedbams
    output:
        set baseName, file("${baseName}.unmapped.marked.bam"),
                      file("${baseName}.unmapped.marked_metrics.tsv") into ch_markedbams

    
    module      'java'
    
    script:
    """
    java -Dpicard.useLegacyParser=false -Xmx30g -jar $picardJar MarkIlluminaAdapters \
        -INPUT $bam \
        -OUTPUT "${baseName}.unmapped.marked.bam" \
        -METRICS "${baseName}.unmapped.marked_metrics.tsv"
    """
}

process alignBwa {
    

    label 'genomics_2'
    
    input:
        set baseName, file(bam), file(metrics) from ch_markedbams
    output:
        set baseName, file("${baseName}.aligned.bam") into ch_mapped, ch_forMetrics1

    publishDir path: './output/intermediate', mode: 'copy'

    module      bwaModule
    module	    'samtools'
    module      'picard'
    
    

    script:
    """
    set -o pipefail
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar SamToFastq \
        -I "$bam" \
        -FASTQ '/dev/stdout' -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 \
        -INTERLEAVE true -NON_PF true -TMP_DIR "$tmp_dir" | \
    bwa mem -M -t ${task.cpus} -p $ref /dev/stdin | \
    java -Dpicard.useLegacyParser=false -Xmx6G -jar $picardJar MergeBamAlignment \
        -ALIGNED_BAM '/dev/stdin' -UNMAPPED_BAM "$bam" \
        -OUTPUT "${baseName}.aligned.bam" -R "$ref" -ADD_MATE_CIGAR true \
        -CLIP_ADAPTERS false -MAX_INSERTIONS_OR_DELETIONS '-1' \
        -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -TMP_DIR "$tmp_dir"
    """
}

process indexPreUmiBam {
        
    label 'genomics_3'
    
    input:
        set baseName, file(bam) from ch_mapped
    output:
        set baseName, file(bam), file("${baseName}.aligned.bam.bai") into ch_indexedMapped
    publishDir path: './output/intermediate', mode: 'copy'

 
    module      'samtools'
    

    script:
    """
    samtools index $bam ${baseName}.aligned.bam.bai
    """

}

ch_tumorPREUMI  = Channel.create()
ch_normalPREUMI = Channel.create()

//split single bam channel into tumor and normal **CURRENTLY RELIES ON "SAMPLE_[FFPE|NORMAL]" naming scheme
//*****FFPE means normal in this study*********
ch_indexedMapped.choice(ch_normalPREUMI, ch_tumorPREUMI){ a -> a[0] =~ /FFPE/ ? 0 : 1 }

//split SAMPLE from FFPE|NORMAL so channels can be joined by sample
ch_normalSplitPREUMI = ch_normalPREUMI.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}
ch_tumorSplitPREUMI = ch_tumorPREUMI.map{ baseName, bam, bai -> [ baseName.split('_')[0], baseName.split('_')[1], bam, bai]}

//merge tumor and normal back together by sample number. 
ch_tumorNormalPairsPREUMI = ch_tumorSplitPREUMI.join(ch_normalSplitPREUMI)


ch_bedSegments2 = Channel.fromPath("$padded_bed").splitText( by: 50000, file: "seg")
ch_vardictPREUMI= ch_tumorNormalPairsPREUMI.combine(ch_bedSegments2)

process runVardictPREUMI {
    
    label 'genomics_3'
    
    input:
        set sample, ttype, file(tbam), file(tbai), ntype, file(nbam), file(nbai), file(segment) from ch_vardictPREUMI
    output:
        set sample, file(tbam), file(nbam), file("${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv") into ch_rawVardictSegmentsPREUMI

      
    script:
    """
    export PATH=/home/jste0021/scripts/git_controlled/vardict_testing/VarDictJava/build/install/VarDict/bin/:$PATH
    VarDict -G ${ref} -f 0.01 -N "${tbam}|${nbam}" \
        -b "${tbam}|${nbam}" -th ${task.cpus} --nosv -c 1 -S 2 -E 3 -g 4 ${segment} \
        > "${sample}.${ttype}_v_${ntype}.${segment}.somatic.vardict.tsv"
    """ 

}


ch_collatedSegmentsPREUMI = ch_rawVardictSegmentsPREUMI.map{ sample, tbam, nbam, segment -> [sample, tbam.name, nbam.name, segment]}.groupTuple(by: [0,1,2])


process catSegmentsPREUMI {
    
    label 'genomics_3'
    
    echo true
    input: 
        set sample, tbam, nbam, file(tsv) from ch_collatedSegmentsPREUMI
    output: 
        set sample, tbam, nbam, file("${sample}.collated.vardict.tsv") into ch_rawVardictPREUMI

    
    
    script:
    
    myfiles = tsv.collect().join(' ')

    """
    cat ${myfiles} > ${sample}.collated.vardict.tsv
    """

}

process makeVCFPREUMI {
    
    label 'genomics_3'
    
    input:
        set sample, tbam, nbam, file(tsv) from ch_rawVardictPREUMI
    output:
        set sample, file("${sample}.somatic.vardict.vcf") into ch_outputVCFPREUMI
    
    //publishDir path: './output/preUMI/intermediate', mode: 'copy'
    
    
    script:
    """  
    module purge
    module load R/3.5.1
    cat $tsv | ome/jste0021/scripts/VarDict-1.7.0/bin/teststrandbias.R | \
        /home/jste0021/scripts/VarDict-1.7.0/bin/var2vcf_valid.pl -N "${tbam}|${nbam}"  \
        -f 0.01 > "${sample}.somatic.vardict.vcf"
    """
}

process reheaderPREUMIVCF {
    
    label 'genomics_3'
    
    input:
        set sample, file(vcf) from ch_outputVCFPREUMI
    output:
        set sample, file("*.vcf.gz") into ch_reheaderVCFPREUMI

    //publishDir path: './output/preUMI/intermediate', mode: 'copy'
    
        
    module      'bcftools/1.8'

    script:
    """
    bcftools annotate -h ~/vh83/reference/genomes/b37/vcf_contig_header_lines.txt -O v ${vcf} | \
        bcftools sort -o ${sample}.somatic.vardict.reheader.vcf.gz -O z -
    """

}

process sortVCFSPREUMI {

    label 'genomics_3'
    
    input:
        set baseName, file(vcf) from ch_reheaderVCFPREUMI
    output:
        set baseName, file("${baseName}.somatic.vardict.reheader.sorted.vcf.gz") into ch_sortedVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'                                    
    
    module      'bcftools/1.8'                                             
    
    script:
    """
    bcftools sort -o "${baseName}.somatic.vardict.reheader.sorted.vcf.gz" -O z ${vcf}
    """
}

process indexVCFSPREUMI {
    
    label 'genomics_3'
    
    input:
        set baseName, file(vcf) from ch_sortedVCFPREUMI
    output:
        set baseName, file(vcf), file("${baseName}.somatic.vardict.reheader.sorted.vcf.gz.tbi") into ch_indexedVCFPREUMI

    publishDir path: './output/preUMI/intermediate', mode: 'copy'                                    
    
    
    module      'bcftools/1.8'                                         
    
    script:
    """
    bcftools index -f --tbi ${vcf} -o ${baseName}.somatic.vardict.reheader.sorted.vcf.gz.tbi
    """
}

process vt_decompose_normalisePREUMI {
        
    label 'genomics_3'
    
    input:
        set baseName, file(vcf), file(tbi) from ch_indexedVCFPREUMI
    output:
        set baseName, file("${baseName}.somatic.vardict.reheader.sorted.vt.vcf.gz") into ch_vtDecomposeVCFPREUMI

    //publishDir path: './output/preUMI/intermediate', mode: 'copy'

    module      'vt'
    

    """
    vt decompose -s $vcf | vt normalize -r $ref -o "${baseName}.somatic.vardict.reheader.sorted.vt.vcf.gz" -
    """
}

process apply_vepPREUMI {

    label 'vep'
    
    input:
        set baseName, file(vcf) from ch_vtDecomposeVCFPREUMI
    output:
        set baseName, file("${baseName}.somatic.vardict.reheader.sorted.vt.vep.vcf") into ch_vepVCFPREUMI

    publishDir path: './output/preUMI/somatic_annotated', mode: 'copy'

        
    module      'vep/90'

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o "${baseName}.somatic.vardict.reheader.sorted.vt.vep.vcf"
    """
}


//ch_forMetrics = ch_forMetrics1.concat(ch_forMetrics2)
//ch_forMetrics1.concat(ch_forMetrics2).into{ch_forMultipleMetrics;ch_forHSMetrics}
ch_forMetrics1.into{ch_forMultipleMetrics;ch_forHSMetrics}

process collectHSMetrics {

    label 'genomics_3'
    
    input:
        set sample, file(bam) from ch_forHSMetrics
    output:
        set sample, file("*.HSmetrics.txt"), file("*.pertarget.txt") into ch_metrics
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    

    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx6G -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${bam.baseName}.HSmetrics.txt" \
        -R ${ref} \
        -BI $panel_int \
        -TI $padded_int \
        --PER_TARGET_COVERAGE "${bam.baseName}.pertarget.txt"
    """
}

process collectMultipleMetrics {

    label 'genomics_3'
    
    input:
        set sample, file(bam) from ch_forMultipleMetrics
    output:
        set sample, file("*multiple_metrics*") into ch_metrics2
    
    publishDir path: './output/metrics/multiple', mode: 'copy'
    
    

    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx6G -jar ${picardJar} CollectMultipleMetrics \
        -I $bam \
        -O ${bam.baseName}.multiple_metrics \
        -R $ref
    """


}
