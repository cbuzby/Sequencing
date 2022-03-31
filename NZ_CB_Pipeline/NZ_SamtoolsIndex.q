#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_index
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

module purge

module load bwa/intel/0.7.17
module load samtools/intel/1.14

samtools index NZ_UnselectedA.sort
