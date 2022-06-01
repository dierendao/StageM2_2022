#!/bin/sh
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=fast
#SBATCH --account=scrna_gmp_rara
#SBATCH --time 24:00:00
#SBATCH --mem=100000

# module loading
module load singularity cellranger/4.0.0

cellranger mkgtf \
  genes_mm10_2020.gtf \
  genes_mm10_2020.filtered.gtf \
  --attribute=gene_type:protein_coding

