#!/bin/sh
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=fast
#SBATCH --account=scrna_gmp_rara
#SBATCH --time 24:00:00
#SBATCH --mem=100000

# add Rscript path in my PATH
export PATH=$PATH:/shared/ifbstor1/software/miniconda/envs/r-4.1.1/bin/

# module loading
#module load singularity cellranger/4.0.0 
module load r/4.1.1

Rscript /shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/R_src/Seurat4_integration.R
