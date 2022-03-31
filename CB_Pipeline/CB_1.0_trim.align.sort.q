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

#java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 \
#	$1 $2.trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#bwa mem $REF HNGLVDRXY_n01_CuSO4_CSSI_$1.fastq > g.$1.aln-se.sam
#bwa mem -Y -K 100000000 -R "@RG\tID:$1\tLB:$1\tPL:ILLUMINA\tPM:NOVASEQ\tSM:$1" \
#       $REF HNGLVDRXY_n01_CuSO4_CSSI_$1.fastq > gR2.$1.aln-se.sam

bwa mem -Y -K 100000000 -R "@RG\tID:$1\tLB:$1\tPL:ILLUMINA\tPM:NOVASEQ\tSM:$1" \
       $REF $1.trimmed.fq.gz > gR3.$1.aln-se.sam

java -jar $PICARD_JAR SortSam \
          INPUT=gR3.$1.aln-se.sam \
          OUTPUT=gR3.$1.aln-se.bam \
          SORT_ORDER=coordinate

