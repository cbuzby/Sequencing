#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_merge_edit

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/
module load bamtools/intel/2.5.1

REF=/scratch/cb4097/Sequencing/Reference/*.fna

#split the files
bamtools split -in $1 -reference

