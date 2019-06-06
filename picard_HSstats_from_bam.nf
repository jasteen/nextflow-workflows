input_path =file("./bams")
panel_int      = file('/projects/vh83/reference/IDT_exome_panel/AML_AMLv3_Targets_Standard.b37.interval_list')
padded_int     = file('/projects/vh83/reference/IDT_exome_panel/AML_AMLv3_Targets_Standard.b37.interval_list')

//Variables
ref=file("/projects/uc23/reference/genomes/bwa_index/human_g1k_v37_decoy.fasta")

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores          = 12
globalMemoryS     = '6 GB'
globalMemoryM     = '32 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '1h'
globalTimeL       = '6h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

Channel
    .fromPath("${input_path}/*.bam")
    .map{ file -> tuple(file.name.take(file.name.lastIndexOf('.')), file) }
    .into { ch_1 }

Channel
    .fromPath("${input_path}/*.bai")
    .map{ file -> tuple(file.name.take(file.name.lastIndexOf('.')), file) }
    .into { ch_2 }


process collectHSMetrics {

    input:
        set sample, file(bam), file(bai) from ch_1.join(ch_2)
    output:
        set sample, file("*.HSmetrics.txt"), file("*.perbase.txt"), file("*.pertarget.txt") into ch_metrics_unused2
    
    publishDir path: './output/metrics/coverage', mode: 'copy'
    
    executor    globalExecutor
    stageInMode globalStageInMode
    cpus        1
    memory      globalMemoryM
    time        globalTimeL
    queue       globalQueueL

    script:

    """
    module purge
    module load R/3.5.1
    java -Dpicard.useLegacyParser=false -Xmx6G -jar ${picardJar} CollectHsMetrics \
        -I ${bam} \
        -O "${sample}.HSmetrics.txt" \
        -R ${ref} \
        -BI $panel_int \
        -TI $padded_int \
        --PER_BASE_COVERAGE "${sample}.perbase.txt" \
        --PER_TARGET_COVERAGE "${sample}.pertarget.txt"
    """
}


