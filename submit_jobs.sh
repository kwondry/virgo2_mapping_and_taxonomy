#!/bin/bash
#SBATCH --job-name=virgo2_mapping_and_taxonomy
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

mkdir -p logs

# Run snakemake with SLURM executor
snakemake \
    --executor slurm \
    --default-resources slurm_account=kwon slurm_partition=short \
    --use-conda \
    --conda-frontend mamba \
    --conda-prefix envs \
    --jobs 500 \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete