#!/usr/bin/env bash

#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=5 --mem=7000
#SBATCH --time=48:00:00
#SBATCH --job-name="chipfreq"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL

PROJECT="exome_aibst"
HOMEDIR="/cbio/users/mamana/exome_aibst/"
OUTDIR="/scratch3/users/mamana/${PROJECT}"

mkdir -p "${OUTDIR}"
cd ${OUTDIR}

nextflow \
    -log ${OUTDIR}/nextflow.log \
    run /cbio/users/mamana/popfreqs/main.nf \
    -c /cbio/users/mamana/exome_aibst/ILIFU/annot_vcf_with_freq.config \
    -resume \
    -profile slurm,singularity