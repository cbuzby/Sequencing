#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=NZ_PipelineAll
#SBATCH --mail-user=cb4097@nyu.edu

cd /scratch/cb4097/Sequencing/HNGLVDRXY/


module purge
module load gatk/4.2.0.0

gatk VariantsToTable \
     -V $1 \
     -F CHROM -F POS -F REF -F ALT \
     -GF AD -GF DP -GF GQ -GF PL \
     -O $1.output.table

