#########################################
# Input data
#########################################
gtex_genotype_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_genotype/"

gtex_expression_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/"

gtex_covariate_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_covariates/"





#########################################
# Output data
#########################################
output_root="/n/scratch3/users/b/bes710/gene_heritability/gtex/"

# Output directory containing estimates of gene h2
gene_h2_estimate_dir=$output_root"gene_h2_estimates/"





##########################################################
# Compute gene h2 for all genes in a single tissue using:
### 1. GREML
### 2. GGVE (gene variance estimation)
###########################################################
cis_window="100000"
tissue_name="Adipose_Subcutaneous"
chrom_num="21"


sh estimate_gene_h2_with_greml_and_ggve_on_single_tissue_chromosome.sh $gtex_genotype_dir $gtex_expression_dir $gtex_covariate_dir $cis_window $tissue_name $chrom_num $gene_h2_estimate_dir
