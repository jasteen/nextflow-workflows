#!/usr/bin/env nextflow

// Required Inputs
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')

//project specific bed files

vardictBed       = file("/projects/vh83/reference/brastrap_specific/vardict/BRASTRAP_721717_8column.bed")
intervalFile     = file("/projects/vh83/reference/brastrap_specific/BRA-STRAP_621717_100.final.roverfile_g37.numsort.sorted.bed")
restrictedBed    = file("/projects/vh83/reference/brastrap_specific/BRA-STRAP_coding_regions_targeted_sort.bed")
primer_bedpe_file= file("/projects/vh83/reference/prostrap/final_prostrap_b37_bedpe_bamclipper.txt")

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

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	      = 4
vepCores          = 4
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '15m'
globalTimeL       = '15m'
globalQueueS      = 'short'
globalQueueL      = 'short'

// Creating channel from input directory
ch_inputFiles = Channel.fromFilePairs("${inputDirectory}/*_R{1,2}_001.fastq.gz")


process align_bwa {

    input:
        set baseName, file(fastqs) from ch_inputFiles
    output:
        set baseName, file("${baseName}.hq.sorted.bam"), file("${baseName}.hq.sorted.bam.bai") into ch_mappedBams

    publishDir path: './bam_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module      samtoolsModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL

    """
    bwa mem -M -t ${task.cpus} -R "@RG\\tID:${baseName}\\tSM:${baseName}\\tPU:lib1\\tPL:Illumina" $ref ${fastqs[0]} ${fastqs[1]}  \
        | samtools view -u -h -q 1 -f 2 -F 4 -F 8 -F 256 - \
        | samtools sort -@ $bwaCores -o "${baseName}.hq.sorted.bam"
    samtools index "${baseName}.hq.sorted.bam" "${baseName}.hq.sorted.bam.bai"
    """
}

ch_mappedBams.into{ch_mappedBam1;ch_mappedBam2;ch_mappedBam3;ch_mappedBam4;ch_mappedBam5;ch_forBamClipper}

process run_bamClipper {
    input:
        set baseName, file(bam), file(bai) from ch_forBamClipper               
    output: 
        set baseName, file("${baseName}.hq.sorted.primerclipped.bam"), file("${baseName}.hq.sorted.primerclipped.bam.bai") into ch_forperBase          
    
    publishDir path: './bamclipper', mode: 'copy'                                    
    
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL 
    module      'gnuparallel/20190122'
    module      'samtools'

    """
    ${bamclipper_exe} -b ${bam} -p ${primer_bedpe_file} -n ${task.cpus}
    """
}

//***magic sample collection right here generate list.txt***
ch_forperBase.into{ch_bamList;ch_bams}

//set one version to a list of filenames of the VCF

ch_bamList.map { it -> it[1].name }
       .collectFile(name: 'list.txt', newLine: true)
       .splitText( by: 10, file: "temp")
       .set {ch_bamList_f}

//set the second to all the files

ch_bams
    .collect()
    .set {ch_all_bams}

//awk 'BEGIN{FS=OFS="\t"}{if($0 ~ /^#/)next;call=0; nocall=0;for(i=10; i<=NF; i++)if($i ~ /^\.\/\.:/)nocall++;else call++;print $1, $2, $4, $5, call, nocall}'

process generatePerbaseMetrics {
    input:
        file list from ch_bamList_f
        file '*' from ch_all_bams               
    output: 
        file("mpileup_${list}.vcf.gz") into ch_mpileupOUT           
    
    publishDir path: './bamclipper', mode: 'copy'                                    
    
    
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryL 
    time        globalTimeM
    queue       globalQueueL 
    module      'samtools'
    cpus        '8'
    module      'bcftools'

    """
    bcftools mpileup --threads ${task.cpus} -Oz -d 250 -B -R ${restrictedBed} -a "FORMAT/DP" -f ${ref} -b ${list} -o mpileup_${list}.vcf.gz

    """
//| bcftools call --threads ${task.cpus} -Oz -m -o mpileup_out.vcf.gz 
}

process sortVCFS {

    input:
        file(vcf) from ch_mpileupOUT
    output:
        file("${vcf}.sorted.vcf.gz") into ch_mpileupsortedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    module      bwaModule
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    script:
    """
    bcftools sort -o "${vcf}.sorted.vcf.gz" -O z ${vcf}
    """
}


process indexVCFS {
    input:
        file(vcf) from ch_mpileupsortedVCF
    output:
        set file(vcf), file("*.tbi") into ch_indexedmpileupVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    module      bwaModule
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    script:
    """
    bcftools index -f --tbi ${vcf} -o ${vcf}.tbi
    """
}

//duplicate ch_indexedVCF
ch_indexedmpileupVCF
    .into{ch_mpileuplist;ch_mpileup_files}
//set one version to a list of filenames of the VCF
ch_mpileuplist.map { it -> it[0].name }
       .collectFile(name: 'list.txt', newLine: true)
       .set {ch_mpileup_list_f}
//set the second to all the files
ch_mpileup_files
    .collect()
    .set {ch_mpileup_all_files}

//feed both to the merge so that the indexes are available to bcftools
process mergeVCFS {
    echo true
    publishDir './variants_merged/', mode: 'copy'
    input:
    file list from ch_mpileup_list_f
    file '*' from ch_mpileup_all_files
    
    output:
    file "merged.mpileup.vcf.gz" into ch_mergedVCF

    module     'bcftools/1.8'
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    memory      globalMemoryM 
    time        globalTimeM
    queue       globalQueueL

    
    script: 
    
    """
    bcftools merge -O z -o "merged.mpileup.vcf.gz" -l list.txt
    """
}



