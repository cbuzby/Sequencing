#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=nf_full
#SBATCH --output=slurm_full.out

module purge
module load nextflow/20.10.0

cd /scratch/cb4097/Sequencing/
nextflow run /scratch/cb4097/Sequencing/NEXTFLOW_TESTING/TrimAlignSortBQSR_All.nf -resume -c /scratch/cb4097/Sequencing/NEXTFLOW_TESTING/nextflow.CSS8_CuSO4_Apr23.config -with-timeline timeline_full.html -with-report report_full.html
