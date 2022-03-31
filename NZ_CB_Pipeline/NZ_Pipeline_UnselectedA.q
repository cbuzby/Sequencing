#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_PipelineAll
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

module purge

module load bwa/intel/0.7.17
module load samtools/intel/1.14
module load bcftools/intel/1.14

REF=/scratch/cb4097/Sequencing/Reference/*.fna

bwa aln -t 12 $REF HNGLVDRXY_n01_CuSO4_CSSI_UnselectedA.fastq > NZ_UnselectedA.sai

bwa samse $REF NZ_UnselectedA.sai HNGLVDRXY_n01_CuSO4_CSSI_UnselectedA.fastq > NZ_UnselectedA.sam

samtools view -bS NZ_UnselectedA.sam | samtools sort -o NZ_UnselectedA.sort -O BAM

samtools index NZ_UnselectedA.sort

bcftools mpileup -f $REF NZ_UnselectedA.sort > NZ__UnselectedA.sort.bcf

bcftools call -c --ploidy 1 NZ__UnselectedA.sort.bcf > NZ_UnselectedA.snps

