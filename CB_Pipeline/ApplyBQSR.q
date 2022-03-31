#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_Pipeline

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17

REF=/scratch/cb4097/Sequencing/Reference/*.fna

gatk ApplyBQSR \
   -R $REF \
   -I $1 \
   --bqsr-recal-file $2 \
   -O ${1}_bqsr.bam
 
