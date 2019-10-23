#!/bin/bash


# The name of the job:
#SBATCH --job-name="heidi"

# Partition for the job:
#SBATCH --partition=comp

# Account for the job
#SBATCH --account=uc23
#SBATCH --qos=normal
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=64000MB

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-00:00:00


#130243-G.bwa.sorted.connor.bam
#130243-G.bwa.sorted.connor.bam.bai
#130243-T.bwa.sorted.connor.bam
#130243-T.bwa.sorted.connor.bam.bai
#130309-G.bwa.sorted.connor.bam
#130309-G.bwa.sorted.connor.bam.bai
#130309-T.bwa.sorted.connor.bam
#130309-T.bwa.sorted.connor.bam.bai
#130313-G.bwa.sorted.connor.bam
#130313-G.bwa.sorted.connor.bam.bai
#130313-T.bwa.sorted.connor.bam
#130313-T.bwa.sorted.connor.bam.bai
#130313-TP.bwa.sorted.connor.bam
#30313-TP.bwa.sorted.connor.bam.bai

ref="/scratch/vh83/projects/cell-free/human_g1k_v37.fasta.gz"
bed="/scratch/vh83/projects/cell-free/cellfree.bed"
export PATH=/projects/uc23/software/VarDict-1.5.8/bin/:$PATH

####

VarDict -G ${ref} -f 0.00001 -N "130243-T|130243-G" \
    -b "130243-T.bwa.sorted.connor.bam|130243-G.bwa.sorted.connor.bam" -th 12 --nosv -c 1 -S 2 -E 3 -g 4 ${bed} \
    > "130243.TvG.somatic.vardict.tsv"

module purge
module load R/3.5.1
cat 130243.TvG.somatic.vardict.tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/testsomatic.R | \
/home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_paired.pl -N "130243-T|130243-G" -f 0.00001 > "${sample}.somatic.vardict.vcf"

#####

VarDict -G ${ref} -f 0.00001 -N "130309-T|130309-G" \
    -b "130309-T.bwa.sorted.connor.bam|130309-G.bwa.sorted.connor.bam" -th 12 --nosv -c 1 -S 2 -E 3 -g 4 ${bed} \
    > "130309.TvG.somatic.vardict.tsv"

module purge
module load R/3.5.1
cat 130309.TvG.somatic.vardict.tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/testsomatic.R | \
/home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_paired.pl -N "130309-T|130309-G" -f 0.00001 > "${sample}.somatic.vardict.vcf"

####

VarDict -G ${ref} -f 0.00001 -N "130313-T|130313-G" \
    -b "130313-T.bwa.sorted.connor.bam|130313-G.bwa.sorted.connor.bam" -th 12 --nosv -c 1 -S 2 -E 3 -g 4 ${bed} \
    > "130313.TvG.somatic.vardict.tsv"

module purge
module load R/3.5.1
cat 130313.TvG.somatic.vardict.tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/testsomatic.R | \
/home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_paired.pl -N "130313-T|130313-G" -f 0.00001 > "${sample}.somatic.vardict.vcf"

######

VarDict -G ${ref} -f 0.00001 -N "130313-TP|130313-G" \
    -b "130313-TP.bwa.sorted.connor.bam|130313-G.bwa.sorted.connor.bam" -th 12 --nosv -c 1 -S 2 -E 3 -g 4 ${bed} \
    > "130313.TPvG.somatic.vardict.tsv"

module purge
module load R/3.5.1
cat 130313.TPvG.somatic.vardict.tsv | /home/jste0021/scripts/VarDict-1.5.8/bin/testsomatic.R | \
/home/jste0021/scripts/VarDict-1.5.8/bin/var2vcf_paired.pl -N "130313-TP|130313-G" -f 0.00001 > "${sample}.somatic.vardict.vcf"

######

