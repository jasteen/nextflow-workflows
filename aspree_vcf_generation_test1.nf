#!/usr/bin/env nextflow

refFolder      = file("/mspm-data/processing/reference")

//project specific bed files

vardictBed       = file("/mspm-data/processing/bed_files/aspree_brastrap_overlap.8col_hg19.bed") // 8 column BED file suitable for vardict
intervalFile     = file("/mspm-data/processing/bed_files/aspree_brastrap_overlap_hg19.bed") //standard bed file of intervals covered


// Getting Reference Files
refBase          = "$refFolder/ucsc.hg19"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")

// random reference stuff
header           = file("/mspm-data/processing/random_files/vcf_contig_header_lines_hg19.txt")
af_thr           = 0.1
rheader          = file("/mspm-data/processing/random_files/Rheader.txt")


ch_bamIN = Channel.fromFilePairs("/home/jste0021/mt-sinai/**/*.{bam,bam.bai}", flat: true)

process run_vardict {

    label 'small_1'

    input:
        set baseName, file(bam), file(bai) from ch_bamIN               
    output: 
        set baseName, file("${baseName}.tsv") into ch_vardictTSV           
    
    script:
    """
    export PATH=/mspm-data/processing/VarDict/bin:$PATH
    VarDict -k 0 -q 15 -G ${ref} -f ${af_thr} -N "${baseName}" -b ${bam} -c 1 -S 2 -E 3 -g 4 ${intervalFile} > "${baseName}.tsv"
    """
}

process makeVCF {
    
    label 'small_1'

    input:
        set baseName, file(tsv) from ch_vardictTSV
    output:
        set baseName, file("${baseName}.vardict.vcf") into ch_vardictVCFs
    
    script:
    """
    cat ${tsv} | /mspm-data/processing/VarDict/bin/teststrandbias.R | \
        /mspm-data/processing/VarDict/bin/var2vcf_valid.pl -N "${baseName}" \
        -f ${af_thr} -E > "${baseName}.vardict.vcf"
    """
}

process reheaderUMIVCF {
    
    label 'small_1'

    input:
        set sample, file(vcf) from ch_vardictVCFs
    output:
        set sample, file("*.vcf.gz") into ch_reheaderVCF

    script:
  
    """
    bcftools annotate -h ${header} -O v ${vcf} | \
        bcftools sort -o ${sample}.vardict.sorted.vcf.gz -O z -
    """
}

process indexVCFS {

    label 'small_1'

    input:
        set sample, file(vcf) from ch_reheaderVCF
    output:
        set file(vcf), file("${sample}.vardict.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir './variants_raw/', mode: 'copy'

    script:
  
    """
    bcftools index -f --tbi ${vcf} -o ${sample}.vardict.sorted.vcf.gz.tbi
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
    file("final_merge.vardict.vcf.gz") into ch_mergedfinalVCF

    script: 
    
    """
    split -l 500 list.txt temp_shorter_list_
    for i in temp_shorter_*; do bcftools merge -m none -l \$i -O z -o \$i.merged.vcf.gz; bcftools index \$i.merged.vcf.gz; done
    ls *merged.vcf.gz > list3.txt
    bcftools merge -R ${intervalFile} -m none -O z -o "final_merge.vardict.vcf.gz" -l list3.txt
    """
}

process vt_decompose_normalise {
        
    label 'small_1'

    input:
        file(vcf) from ch_mergedfinalVCF
    output:
        file("vardict.merged.vt.vcf.gz") into ch_vtDecomposeVCF
    
    script:
    """
    vt decompose -s $vcf | vt normalize -r $ref -o "vardict.merged.vt.vcf.gz" -
    """
}

