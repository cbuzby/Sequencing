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

bcftool mpileup -Ou -f $REF NZ_UnselectedA.sort | bcftools call -cv --ploidy 1 -Ob -o NZ_UnselectedA.calls.bcf
