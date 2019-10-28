#!/usr/bin/env nextflow

refFolder      = file("~/reference")

//project specific bed files

vardictBed       = file("") // 8 column BED file suitable for vardict
intervalFile     = file("~/WG_IAD49736_v2.20131226.designed.bed") //standard bed file of intervals covered


// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")


millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")

// random reference stuff
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")


ch_bamIN = Channel.fromFilePairs("./mt-sinai/**/*.{bam,bam.bai}")

process run_vardict {

    input:
        set baseName, file(bam), file(bai) from ch_mappedBam1               
    output: 
        set baseName, file("${baseName}.tsv") into ch_vardictTSV           
    
    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    """
    export PATH=/home/jste0021/scripts/VarDict-1.5.8/bin/:$PATH
    VarDict -G ${ref} -f 0.1 -N "${baseName}" -b ${bam} -c 1 -S 2 -E 3 -g 4 ${vardictBed} > "${baseName}.tsv"
    """
}

process makeVCF {
    input:
        set baseName, file(tsv) from ch_vardictTSV
    output:
        set baseName, file("${baseName}.vardict.vcf") into ch_vardictVCFs
    
    script:

    """
    module purge
    module load R/3.5.1
    cat ${tsv} | /home/jste0021/scripts/VarDict-1.5.8/bin/teststrandbias.R | \
        /home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_valid.pl -N "${baseName}" \
        -f 0.1 -E > "${baseName}.vardict.vcf"
    """
}

process reheaderUMIVCF {
    
    label 'small_1'

    input:
        set sample, file(vcf) from ch_vardictVCFS
    output:
        set sample, file("*.vcf.gz") into ch_reheaderVCF

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
        set sample, file(vcf), file("${sample}.vardict.sorted.vcf.gz.tbi") into ch_indexedVCF

    module     'bcftools/1.8'

    script:
  
    """
    bcftools index -f --tbi ${vcf} -o ${sample}.UMI.reheader.sorted.vcf.gz.tbi
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

    publishDir path: './output/vcf/', mode: 'copy', pattern: "*.vcf"
    
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
                      -o "${sample}.UMI.reheader.sorted.vt.vep.vcf"
    """
}
