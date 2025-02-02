#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_samse
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

module purge

module load bwa/intel/0.7.17

REF=/scratch/cb4097/Sequencing/Reference/*.fna

bwa samse $REF NZ_UnselectedA.sai HNGLVDRXY_n01_CuSO4_CSSI_UnselectedA.fastq > NZ_UnselectedA.sam

