#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=couldrun


module purge
module load bcftools/intel/1.14
module load gatk/4.2.0.0


for i in ${1}*vcf; do bgzip -c $i > ${i}.gz; done
echo ${1}*vcf.gz |  xargs -n1 tabix -p vcf
bcftools concat -o unsortedcat.vcf -a -D ${1}*vcf.gz
bcftools sort -Oz -o ${1}.SortedCat.vcf unsortedcat.vcf

myfile=${1}.SortedCat.vcf

gatk VariantsToTable \
     -V ${myfile} \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O ${myfile}.output.table

