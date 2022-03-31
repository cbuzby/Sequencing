#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_merge_edit

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/

REF=/scratch/cb4097/Sequencing/Reference/*fna

#$1 is a bam to be called

echo $REF
echo $1

#gatk HaplotypeCaller -I $1 -R $REF -ploidy 1 -O ${1}.g.vcf.gz -ERC GVCF

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I $1 \
   -O ${1}.g.vcf.gz \
   -ERC GVCF

#gatk --java-options "-Xmx4g" GenotypeGVCFs \
#   -R  $REF \
#   -V ${1}.g.vcf.gz \
#   -O ${1}.GVCF_output.vcf.gz

