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

bwa aln -t 12 $REF *$1*.fastq > NZ.$1.sai

bwa samse $REF NZ.$1.sai *$1*.fastq > NZ.$1.sam

samtools view -bS NZ.$1.sam | samtools sort -o NZ.$1.sort -O BAM

samtools index NZ.$1.sort

bcftools mpileup -f $REF NZ.$1.sort > NZ.$1.sort.bcf

bcftools call -c --ploidy 1 NZ.$1.sort.bcf > NZ.$1.snps
#adding -f GQ to call for purposes of getting GQ
