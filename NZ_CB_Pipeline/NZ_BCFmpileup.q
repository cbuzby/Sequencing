#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_bcfmpileup
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

module purge

module load bwa/intel/0.7.17
module load samtools/intel/1.14
module load bcftools/intel/1.14

REF=/scratch/cb4097/Sequencing/Reference/*.fna

#samtools mpileup -I -uf $REF NZ_UnselectedA.sort > NZ_UnselectedA.sort.bcf

bcftools mpileup -f $REF NZ_UnselectedA.sort > NZ__UnselectedA.sort.bcf
