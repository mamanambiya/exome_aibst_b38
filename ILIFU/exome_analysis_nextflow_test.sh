#!/usr/bin/env bash

#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=2 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="exome"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
#SBTACH -o /scratch/users/mamana/exome_aibst/LOG/.nextflow.out

cd /cbio/users/mamana/exome_aibst_test

# nextflow ~/exome_aibst/main.nf -c ~/exome_aibst/ILIFU/exome_analysis_nextflow_test.config -resume -profile singularity
nextflow ~/exome_aibst/main.nf -c ~/exome_aibst/ILIFU/exome_analysis_nextflow_test.config -resume -profile singularity,slurm
