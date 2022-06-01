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
cellranger mkref --genome=REFERENCE_plasmide \
        --fasta=mm10_plasmide_PLZF_RARA.fa \
        --genes=mm10_plasmide_PLZF_RARA.gtf

