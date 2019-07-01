#!/usr/bin/env nextflow

// Required Inputs

project_name   = "PALB2_isoform_analysis"
refFolder      = file("/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/")
inputDirectory = file('./fastqs')
tmp_dir        = file('/scratch/vh83/tmp/')
chromsizes     = file('/scratch/vh83/projects/small_projects/palb2_nanopore/chrom_sizes.txt')
humanGTF       = file('/scratch/vh83/projects/small_projects/palb2_nanopore/Homo_sapiens.GRCh37.87.chr.gtf')

// Getting Reference Files
refBase          = "$refFolder/human_g1k_v37_decoy"
ref              = file("${refBase}.fasta")
refDict          = file("${refBase}.dict")
refFai           = file("${refBase}.fasta.fai")

// Tools
condaModule            = 'miniconda3/4.1.11-python3.5'          
samtoolsModule         = 'samtools/1.9-gcc5'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	      = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

ch_fastqs = Channel.fromPath("${inputDirectory}/*.fastq").collect()

process catFastqs {
    
    publishDir path: './output', mode: 'copy'
    
    input:
        file '*.fastq' from ch_fastqs
    output:
        file "${project_name}.all.fastq" into ch_catFastq, ch_catFastq2

    publishDir path: './output', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        '1h'
    queue       globalQueueL

    script:
    """
    cat *.fastq > "${project_name}.all.fastq"
    """
}


process flairAlign {
    
    publishDir path: './output', mode: 'copy'
    
    input:
        file(fastq) from ch_catFastq
    output:
        set file("${project_name}.flair_aligned.sam"), file("${project_name}.flair_aligned.bed")  into ch_align, ch_align2

    publishDir path: './output', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        8
    module      condaModule
    conda       '/home/jste0021/.conda/envs/py3.5/'
    module      samtoolsModule
    memory      globalMemoryM
    time        '6h'
    queue       globalQueueL

    script:
    """
    python ~/scripts/git_controlled/flair/flair.py align -v1.3 -g ${ref} -r ${fastq} -t ${task.cpus} \
           -o "${project_name}.flair_aligned"
    """
}

process makeBam{

    publishDir path: './output', mode: 'copy'
    
    input:
        set file(sam), file(bed) from ch_align2
    output:
        set file("*.bam"), file("*.bai") into ch_bams

    publishDir path: './output', mode: 'copy'
    
    cache         'lenient'
    executor      globalExecutor
    stageInMode   globalStageInMode
    cpus          1
    module        samtoolsModule
    memory        globalMemoryM
    time          '6h'
    queue         globalQueueL
    errorStrategy 'ignore'

    script:
    """
    samtools view -bS ${sam} | sort -o "${sam}.bam"
    samtools index "${sam}.bam"
    """
}


process flairCorrect {
    
    publishDir path: './output', mode: 'copy'
    
    input:
        set file(sam), file(bed) from ch_align
    output:
        set file("*_all_corrected.bed"), file("*_all_inconsistent.bed"), file("*.psl") into ch_correct

    publishDir path: './output', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        8
    module      condaModule
    conda       '/home/jste0021/.conda/envs/py3.5/'
    module      samtoolsModule
    memory      globalMemoryM
    time        '6h'
    queue       globalQueueL

    script:
    """
    python ~/scripts/git_controlled/flair/flair.py correct -g ${ref} -q ${bed} -f $humanGTF -c $chromsizes -t ${task.cpus}
    """
}


process flairCollapse {
    
    publishDir path: './output', mode: 'copy'
    
    input:
        set file(corrected), file(inconsistent), file(psl) from ch_correct
        file(fastq) from ch_catFastq2
    output:
        set file("*isoforms.psl"), file("*isoforms.gtf"), file("*isoforms.fa") into ch_collapse

    publishDir path: './output', mode: 'copy'
    
    cache       'lenient'
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        8
    module      condaModule
    conda       '/home/jste0021/.conda/envs/py3.5/'
    module      samtoolsModule
    memory      globalMemoryM
    time        '6h'
    queue       globalQueueL

    script:
    """
    python ~/scripts/git_controlled/flair/flair.py collapse --keep_intermediate -g ${ref} -r ${fastq} -q ${psl} -t ${task.cpus}
    """
}