#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_merge_edit

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/
module load bamtools/intel/2.5.1

REF=/scratch/cb4097/Sequencing/Reference/*.fna

echo $REF

echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}
echo "HKTFTDRX2 is the new merged file name"

#samtools merge HNGLVDRXY_gR_merged.bam gR*.bam
#samtools merge HNGLVDRXY_${1}_merged.bam $2

samtools merge HVYTYDRX2 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}
samtools index HVYTYDRX2

#split the files
bamtools split -in HVYTYDRX2 -reference

#Don't run this; do it in parallel because of timing:
#gatk HaplotypeCaller -I $2 -R $REF -ploidy 1 -O HNGLVDRXY_${1}_merged.vcf
