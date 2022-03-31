#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_Pipeline

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17

REF=/scratch/cb4097/Sequencing/Reference/*.fna

gatk BaseRecalibrator -R $REF -I $1 --known-sites $2 -O ${1}_recal_data.table

##Old Version
#java -jar $GATK_JAR -T BaseRecalibrator -R GCF_000001405.33_GRCh38.p7_chr20_genomic.fna -I dedup_reads.bam -knownSites 1000G_omni2.5.hg38.vcf.gz -o recal_data.table
