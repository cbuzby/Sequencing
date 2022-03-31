#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_viewsort
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

module purge

module load bwa/intel/0.7.17
module load samtools/intel/1.14

REF=/scratch/cb4097/Sequencing/Reference/*.fna

samtools view -bS NZ_UnselectedA.sam | samtools sort -o NZ_UnselectedA.sort -O BAM
#samtools view -bS NZ_UnselectedA.sam | samtools sort - NZ_Unselected.sort

