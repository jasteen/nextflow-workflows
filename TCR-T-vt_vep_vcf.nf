refFolder        = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")



//Annotation resources
dbsnp_b37       = file("/projects/vh83/reference/genomes/b37/accessory_files/dbsnp_138.b37.vcf")
other_vep       = file("/usr/local/vep/90/ensembl-vep/cache")
vep_brcaex      = file("/projects/vh83/reference/annotation_databases/BRCA-Exchange/BRCA-exchange_accessed-180118/BRCA-exchange_accessed-180118.sort.vcf.gz")
vep_gnomad      = file("/projects/vh83/reference/annotation_databases/gnomAD/gnomad.exomes.r2.0.2.sites/gnomad.exomes.r2.0.2.sites.vcf.gz")
vep_revel       = file("/projects/vh83/reference/annotation_databases/REVEL/REVEL-030616/revel_all_chromosomes.vcf.gz")
vep_maxentscan  = file("/projects/vh83/reference/annotation_databases/MaxEntScan/MaxEntScan_accessed-240118")
vep_exac        = file("/projects/vh83/reference/annotation_databases/ExAC/ExAC_nonTCGA.r0.3.1/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz")
vep_dbnsfp      = file("/projects/vh83/reference/annotation_databases/dbNSFP/dbNSFPv2.9.3-VEP/dbNSFP-2.9.3.gz")
vep_dbscsnv     = file("/projects/vh83/reference/annotation_databases/dbscSNV/dbscSNV1.0-VEP/dbscSNV.txt.gz")
vep_cadd        = file("/projects/vh83/reference/annotation_databases/CADD/CADD-v1.3/1000G_phase3.tsv.gz")

// Global Resource Configuration Options

ch_vcfIN = Channel.fromPath("*.vcf.gz").map{file -> tuple(file.name.take(file.name.lastIndexOf('.')), file)}


process indexVCFS {
  
    label 'small_1'

    input:
        set baseName, file(vcf) from ch_vcfIN
    output:
        set baseName, file(vcf), file("${baseName}.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'
    
    script:
    """
    bcftools index -f --tbi ${vcf} -o ${baseName}.vcf.gz.tbi
    """
}

process vt_decompose_normalise {

    label 'medium_6h'
        
    input:
        set baseName, file(vcf), file(tbi) from ch_indexedVCF
    output:
        set baseName, file("${baseName}.vt.vcf.gz") into ch_vtDecomposeVCF

    publishDir path: './variants_raw_out', mode: 'copy'

    
    module      'vt'
    
    """
    vt decompose -s $vcf | vt normalize -n -r $ref -o "${baseName}.vt.vcf.gz" -
    """
}

process apply_vep {

    label 'vep'

    input:
        set baseName, file(vcf) from ch_vtDecomposeVCF
    output:
        set baseName, file("${baseName}.vt.vep.vcf") into ch_vepVCF

    publishDir path: './annotated_variants', mode: 'copy'

    
    
    module      'vep/90'

    """
    vep --cache --dir_cache $other_vep \
                      --assembly GRCh37 --refseq --offline \
                      --fasta $ref \
                      --sift b --polyphen b --symbol --numbers --biotype \
                      --total_length --hgvs --format vcf \
                      --vcf --force_overwrite --flag_pick --no_stats \
                      --uniprot --protein --ccds --canonical --domains \
                      --fork ${task.cpus} \
                      -i ${vcf} \
                      -o "${baseName}.vt.vep.vcf"
    """
}

