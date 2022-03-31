#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_merge_edit

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/

REF=/scratch/cb4097/Sequencing/Reference/*.fna

echo $REF

gatk ApplyBQSR \
   -R $REF \
   -I ${1}.bam \
   --bqsr-recal-file ${2}recalibration.table \
   -O ${1}.bqsr.output.bam
