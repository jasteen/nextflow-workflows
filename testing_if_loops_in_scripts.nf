#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./bams')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files

params.vardictBed       = ""
params.intervalFile     = ""
params.restrictedBed    = ""

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")
millsIndels      = file("${refFolder}/accessory_files/Mills_and_1000G_gold_standard.indels.b37.vcf")
dbSNP            = file("${refFolder}/accessory_files/dbsnp_138.b37.vcf")
genome_file      = file("/projects/vh83/reference/genomes/b37/accessory_files/human_g1k_v37_decoy_GenomeFile.txt")
header           = file("/home/jste0021/vh83/reference/genomes/b37/vcf_contig_header_lines.txt")
af_thr           = 0.1
rheader          = file("/projects/vh83/pipelines/code/Rheader.txt")

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

// Tools
picardJar      = '~/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'
bamclipper_exe = '/projects/vh83/local_software/bamclipper/bamclipper.sh'


// Creating channel from input directory
ch_inputFiles = Channel.fromPath("${inputDirectory}/*.bam").map{file -> tuple(file.name.take(file.name.lastIndexOf('.')), file)}

process generate_bam_index {
    label 'small_short'

    input:
        set baseName, file(bam) from ch_inputFiles             
    output: 
        set baseName, file(bam), file("${baseName}.bam.bai") into ch_mappedBams           
    
    module samtoolsModule

    script:
    """
    samtools index $bam "${baseName}.bam.bai"
    """

}

ch_mappedBams.into{ch_mappedBam1;ch_mappedBam2;ch_mappedBam3;ch_mappedBam4;ch_mappedBam5;ch_mappedBam6}


ch_mappedBam1.branch {
        halo: it[0].contains('locatit')
        hiplex: true
    }
    .set { ch_mappedBam_split }

ch_mappedBam_split.halo.view{"this should be a halo line: $it[0]"}
ch_mappedBam_split.hiplex.view{"this should be a hiplex line: $it[0]"}


/*
process run_vardict_halo {

    label 'vardict_small'

    input:
        set baseName, file(bam), file(bai) from ch_mappedBam1               
    output: 
        set baseName, file("${baseName}.tsv") into ch_vardict_halo_TSV           
    
    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    script:
    """
    export PATH=/home/jste0021/scripts/git_controlled/vardict_testing/VarDictJava/build/install/VarDict/bin/:$PATH
    VarDict -G ${ref} -f 0.1 -N "${baseName}" -p --nosv -b ${bam} -c 1 -S 2 -E 3 -g 4 ${params.intervalFile} > "${baseName}.tsv"
    """
}

process run_vardict_hiplex {

    label 'vardict_small'

    input:
        set baseName, file(bam), file(bai) from ch_mappedBam1               
    output: 
        set baseName, file("${baseName}.tsv") into ch_vardict_hiplex_TSV           
    
    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    script:
    """
    export PATH=/home/jste0021/scripts/git_controlled/vardict_testing/VarDictJava/build/install/VarDict/bin/:$PATH
    VarDict -G ${ref} -f 0.1 -N "${baseName}" -p --nosv -b ${bam} -c 1 -S 2 -E 3 -g 4 ${params.vardictBed} > "${baseName}.tsv"
    """
}

ch_all_TSV = ch_vardict_halo_TSV.join(ch_vardict_hiplex_TSV).view{}

*/

