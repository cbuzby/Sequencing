#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_EP_mergebackvcfs
module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/

cd /scratch/cb4097/Sequencing/HNGLVDRXY_Feb/

samtools index $1
