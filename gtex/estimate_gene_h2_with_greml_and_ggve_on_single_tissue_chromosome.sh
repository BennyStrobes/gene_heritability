#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2GB                         # Memory total in MiB (for all cores)



gtex_genotype_dir="$1"
gtex_expression_dir="$2"
gtex_covariate_dir="$3"
cis_window="$4"
tissue_name="$5"
chrom_num="$6"
gene_h2_estimate_dir="$7"

if false; then
source ~/.bash_profile
fi


python3 estimate_gene_h2_with_greml_and_ggve_on_single_tissue_chromosome.py $gtex_genotype_dir $gtex_expression_dir $gtex_covariate_dir $cis_window $tissue_name $chrom_num $gene_h2_estimate_dir
