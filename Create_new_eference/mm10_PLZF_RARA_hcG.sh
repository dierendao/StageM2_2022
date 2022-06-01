#!/bin/sh
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=fast
#SBATCH --account=scrna_gmp_rara
#SBATCH --time 24:00:00
#SBATCH --mem=100000

# module loading
module load singularity cellranger/4.0.0



# make reference
cellranger mkref --genome=REFERENCE_hcG_utr_3 \
        --fasta=mm10_PLZF_RARA_hcG_complete.fa \
        --genes=mm10_PLZF_RARA_hcG_complete.gtf

