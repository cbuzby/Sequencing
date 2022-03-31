#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_Pipeline

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17

REF=/scratch/cb4097/Sequencing/Reference/*.fna

#Trim data; input is file, then name
#java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 \
#	$1 $2.trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Align to reference genome
bwa mem -Y -K 100000000 -R "@RG\tID:$1\tLB:$1\tPL:ILLUMINA\tPM:NOVASEQ\tSM:$1" \
       $REF ../HNGLVDRXY/$1.trimmed.fq.gz > $1.aln-se.sam

#Sort files
java -jar $PICARD_JAR SortSam \
          INPUT=$1.aln-se.sam \
          OUTPUT=$1.aln-se.bam \
          SORT_ORDER=coordinate

#2/8/21 Update: want to include the BGsomething filter for quality
#also want to not have to re-trim every time, and want to just use the names of those
#also also want to go to a different folder from where the trimmed data is?
gatk BQSR $1.aln-se.bam 
