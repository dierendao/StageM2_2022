#!/bin/sh
#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=fast
#SBATCH --account=scRNA_HSPC_Aging
#SBATCH --time 24:00:00
#SBATCH --mem=100000

# module loading
module load singularity cellranger/4.0.0

# Gene Expression Analysis
cellranger count --id=count_mm10_PLZF_RARA_plasmid_utr_3 \
		   --transcriptome=/shared/ifbstor1/projects/scrna_gmp_rara/stage_babacar/scranseq/reference/CONTROL/REFERENCE_RARA/REFERENCE_plasmide  \
                   --fastqs=/shared/projects/scrna_gmp_rara/GMP_PLZF_RARA/outs/fastq_path  \
                   --sample=Mathilde-23_sCell_cDNA-Ctrl  \
                   --localcores=28 
                 
                   


