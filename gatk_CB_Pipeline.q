#!/bin/bash
#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_Pipeline

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17

cd /scratch/cb4097/Sequencing/HNGLVDRXY/

REF=/scratch/cb4097/Sequencing/Reference/*.fna

#bwa mem $REF HNGLVDRXY_n01_CuSO4_CSSI_$1.fastq > g.$1.aln-se.sam

#gatk MarkDuplicatesSpark -I g.$1.aln-se.sam -M g.$1.dedup_metrics.txt -O g.$1.sorted_dedup_reads.bam

gatk HaplotypeCaller -R $REF -I g.$1.sorted_dedup_reads.bam -o g.$1.variants.vcf

#track changes in github