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

#bcftools mpileup -f $REF -a FORMAT/AD NZFirst/NZ.$1.sort > NZ.$1.sort.bcf
bcftools mpileup -f $REF -a INFO/AD NZFirst/NZ.$1.sort > NZ.$1.mpileup.vcf
bcftools call -m -a gq -v --ploidy 1 NZFirst/NZ.$1.sort.bcf > NZ.$1.call.vcf

#Added Later, optional
bcftools norm -Ou -f $REF -d all -o NZ.$1.snps.norm NZ.$1.snps.vcf
bcftools filter -Ob -e 'QUAL<40 || DP<10' -o NZ.$1.snps.filter 
bcftools index NZ.$1.snps.filter
