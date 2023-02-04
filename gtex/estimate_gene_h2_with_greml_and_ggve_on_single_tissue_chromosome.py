import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import pandas as pd
import pyreadr
import gzip
import time
from sklearn.linear_model import LinearRegression



def load_in_genotype_data(gtex_genotype_dir, tissue_name, chrom_num):
	genotype_stem = gtex_genotype_dir + tissue_name + '_GTEx_v8_genotype_EUR_' + chrom_num
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Centimorgan distances
	G_obj_cm = np.asarray(G_obj.cm)
	# Sample ids
	G_obj_samp_id = np.asarray(G_obj.sample)

	genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm, 'sample_id': G_obj_samp_id}


	return genotype_obj


def load_in_expression_data(gtex_expression_dir, tissue_name, chrom_num):
	expr_mat_file = gtex_expression_dir + tissue_name + '_normalized_expression_matrix_eqtl_ready_chr' + chrom_num + '.txt'
	expr_gene_loc_file = gtex_expression_dir + tissue_name + '_gene_location_matrix_eqtl_ready_chr' + chrom_num + '.txt'

	expr_mat_raw = np.loadtxt(expr_mat_file, dtype=str, delimiter='\t')
	sample_ids = expr_mat_raw[0,1:]
	gene_ids = expr_mat_raw[1:,0]
	expr_mat = expr_mat_raw[1:,1:]

	gene_pos_raw = np.loadtxt(expr_gene_loc_file, dtype=str, delimiter='\t')

	if np.array_equal(gene_pos_raw[1:,0], gene_ids) == False:
		print('assumption eroror')
		pdb.set_trace()

	return expr_mat.astype(float), sample_ids, gene_ids, gene_pos_raw[1:,2].astype(int)


def load_in_covariate_data(gtex_covariate_dir, tissue_name):
	cov_file = gtex_covariate_dir + tissue_name + '_covariates.txt'
	cov_data_raw = np.loadtxt(cov_file, dtype=str, delimiter='\t')
	cov_names = cov_data_raw[1:,0]
	ind_names = cov_data_raw[0,1:]
	cov_data = cov_data_raw[1:,1:].astype(float)
	return cov_data, cov_names, ind_names

def extract_h2_from_greml_hsq_file(hsq_file):
	f = open(hsq_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if line.startswith('V(G)/Vp'):
			hsq = float(data[1])
			hsq_se = float(data[2])
		if line.startswith('Pval'):
			pval = float(data[1])
	f.close()
	return hsq, hsq_se, pval

def run_greml_no_covariate_h2_analysis(chrom_num, gene_tss, expr_vec, genotype_stem, cis_window, expr_ind_ids, tmp_output_stem):
	gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
	# Create gene pheno file
	gene_pheno_file = tmp_output_stem + 'gene_pheno'
	n_samples = len(expr_vec)
	gene_pheno_data = np.hstack((np.zeros((n_samples,1)).astype(int).astype(str), expr_ind_ids.reshape(n_samples,1), expr_vec.astype(str).reshape(n_samples,1)))
	np.savetxt(gene_pheno_file, gene_pheno_data, fmt="%s", delimiter='\t')

	# Run PLINK to get plink file specifically consisting of cis snps
	start_pos = gene_tss - int(cis_window)
	end_pos = gene_tss + int(cis_window)
	plink_window_stem = tmp_output_stem + 'window_plink'
	command_string = 'plink --bfile ' + genotype_stem + ' --keep-allele-order --pheno ' + gene_pheno_file + ' --make-bed --out ' + plink_window_stem + ' --keep ' + gene_pheno_file + ' --chr ' + chrom_num + ' --from-bp ' + str(start_pos) + ' --to-bp ' + str(end_pos) +' --allow-no-sex'
	os.system(command_string)

	# MAKE GRM WITH PLINK
	command_string = 'plink --allow-no-sex --bfile ' + plink_window_stem + ' --make-grm-bin --out ' + plink_window_stem
	os.system(command_string)

	# estimate heritability with GREML
	greml_h2_res_file = tmp_output_stem + 'h2_res' 
	#arg = paste( gcta_path ," --grm ",temp_tissue_specific_stem," --pheno ",raw.pheno.file," --qcovar ",covariate_file," --out ",temp_tissue_specific_stem," --reml --reml-no-constrain --reml-lrt 1",sep='')
	command_string = gcta_path + ' --grm ' + plink_window_stem + ' --pheno ' + gene_pheno_file + ' --out ' + greml_h2_res_file + ' --reml --reml-no-constrain --reml-lrt 1'
	os.system(command_string)

	# Now extract heritabilities
	hsq, hsq_se, hsq_p = extract_h2_from_greml_hsq_file(greml_h2_res_file + '.hsq')

	# Clear temporary files from directory
	os.system('rm ' + tmp_output_stem + '*')

	return hsq, hsq_se, hsq_p

def run_greml_h2_analysis(chrom_num, gene_tss, expr_vec, genotype_stem, cov_mat, cis_window, expr_ind_ids, tmp_output_stem):
	gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
	# Create gene pheno file
	gene_pheno_file = tmp_output_stem + 'gene_pheno'
	n_samples = len(expr_vec)
	gene_pheno_data = np.hstack((np.zeros((n_samples,1)).astype(int).astype(str), expr_ind_ids.reshape(n_samples,1), expr_vec.astype(str).reshape(n_samples,1)))
	np.savetxt(gene_pheno_file, gene_pheno_data, fmt="%s", delimiter='\t')

	# Create covariate file
	gene_covar_file = tmp_output_stem + 'gene_covar'
	n_samples = len(expr_vec)
	gene_covar_data = np.hstack((np.zeros((n_samples,1)).astype(int).astype(str), expr_ind_ids.reshape(n_samples,1), np.transpose(cov_mat).astype(str)))
	np.savetxt(gene_covar_file, gene_covar_data, fmt="%s", delimiter='\t')

	# Run PLINK to get plink file specifically consisting of cis snps
	start_pos = gene_tss - int(cis_window)
	end_pos = gene_tss + int(cis_window)
	plink_window_stem = tmp_output_stem + 'window_plink'
	command_string = 'plink --bfile ' + genotype_stem + ' --keep-allele-order --pheno ' + gene_pheno_file + ' --make-bed --out ' + plink_window_stem + ' --keep ' + gene_pheno_file + ' --chr ' + chrom_num + ' --from-bp ' + str(start_pos) + ' --to-bp ' + str(end_pos) +' --allow-no-sex'
	os.system(command_string)

	# MAKE GRM WITH PLINK
	command_string = 'plink --allow-no-sex --bfile ' + plink_window_stem + ' --make-grm-bin --out ' + plink_window_stem
	os.system(command_string)

	# estimate heritability with GREML
	greml_h2_res_file = tmp_output_stem + 'h2_res' 
	#arg = paste( gcta_path ," --grm ",temp_tissue_specific_stem," --pheno ",raw.pheno.file," --qcovar ",covariate_file," --out ",temp_tissue_specific_stem," --reml --reml-no-constrain --reml-lrt 1",sep='')
	command_string = gcta_path + ' --grm ' + plink_window_stem + ' --pheno ' + gene_pheno_file + ' --qcovar ' + gene_covar_file + ' --out ' + greml_h2_res_file + ' --reml --reml-no-constrain --reml-lrt 1'
	os.system(command_string)

	# Now extract heritabilities
	hsq, hsq_se, hsq_p = extract_h2_from_greml_hsq_file(greml_h2_res_file + '.hsq')

	# Clear temporary files from directory
	os.system('rm ' + tmp_output_stem + '*')

	return hsq, hsq_se, hsq_p

def residualize_gene_expression(expr_vec, cov_mat):
	reg = LinearRegression().fit(np.transpose(cov_mat), expr_vec)
	resid_expr_vec = expr_vec - reg.predict(np.transpose(cov_mat))
	standardized_resid_expr_vec = (resid_expr_vec - np.mean(resid_expr_vec))/np.std(resid_expr_vec)
	return standardized_resid_expr_vec

def standardize_columns_of_matrix(X):
	ncol = X.shape[1]
	X_stand = np.copy(X)

	for col_iter in range(ncol):
		X_stand[:, col_iter] = (X[:, col_iter] - np.mean(X[:, col_iter]))/np.std(X[:, col_iter])
	return X_stand

def multivariate_beta_update(X, y, tau, residual_prec):
	S = np.linalg.inv(np.diag(np.ones(X.shape[0])*tau) + residual_prec*np.dot(X, np.transpose(X)))
	mu = residual_prec*np.dot(np.dot(S, X), y)
	return mu, S

def full_multivariate_beta_update(X, y, tau, residual_prec):
	S = np.linalg.pinv((tau + residual_prec)*np.dot(X, np.transpose(X)))
	mu = residual_prec*np.dot(np.dot(S, X), y)
	return mu, S

# Multivariate learn tau and beta
def full_inference_fixed_residual_var(X_centered, y, residual_variance, max_iter=200):
	expected_tau = 1e-5
	num_snps = X_centered.shape[0]
	num_indi = X_centered.shape[1]
	# Precompute X_X_t
	X_X_t = np.dot(X_centered, np.transpose(X_centered))
	# Iterate
	for global_iter in range(max_iter):
		# Update effects of X on Y
		multivariate_mu, cov = full_multivariate_beta_update(X_centered, y, expected_tau, 1.0/residual_variance)
		# Update Tau
		tau_a = num_indi/2.0
		squared_expectations = cov + np.dot(multivariate_mu.reshape(num_snps,1), multivariate_mu.reshape(1,num_snps))
		tau_b = 0.5*np.trace(np.dot(squared_expectations, X_X_t))
		expected_tau = tau_a/tau_b
		#print(1.0/expected_tau)

	return 1.0/expected_tau

# Multivariate learn tau and beta
def full_inference(X_centered, y, max_iter=200):
	# Initialize some quantities
	expected_tau = 1e-5
	residual_variance = 1.0
	# Precompute some quantities
	num_snps = X_centered.shape[0]
	num_indi = X_centered.shape[1]
	# Precompute X_X_t
	X_X_t = np.dot(X_centered, np.transpose(X_centered))
	
	# Iterate VB
	for global_iter in range(max_iter):
		# Update effects of X on Y
		multivariate_mu, cov = full_multivariate_beta_update(X_centered, y, expected_tau, 1.0/residual_variance)
		# Update Tau
		tau_a = num_indi/2.0
		squared_expectations = cov + np.dot(multivariate_mu.reshape(num_snps,1), multivariate_mu.reshape(1,num_snps))
		tau_b = 0.5*np.trace(np.dot(squared_expectations, X_X_t))
		if tau_b == 0:
			expected_tau_var = 0.0
			break
		expected_tau = tau_a/tau_b
		# Update residual variance
		resid_a = num_indi/2.0
		resid_b = (np.sum(np.square(y) - 2.0*y*np.dot(multivariate_mu, X_centered)) + np.trace(np.dot(squared_expectations, X_X_t)))/2.0

		residual_variance = resid_b/resid_a
		expected_tau_var = 1.0/expected_tau


	return expected_tau_var, residual_variance


def remove_high_ld_snps(X, thresh):
	abs_corr = np.square(np.corrcoef(np.transpose(X)) - np.eye(X.shape[1]))
	n_snps = abs_corr.shape[0]
	valid_snps = []
	for snp_iter in range(n_snps):
		#corry = abs_corr[:snp_iter, snp_iter]
		if len(valid_snps) == 0:
			valid_snps.append(snp_iter)
		else:
			if np.max(abs_corr[snp_iter,valid_snps]) < thresh:
				valid_snps.append(snp_iter)
	valid_snps = np.asarray(valid_snps)
	return X[:,valid_snps]

def run_ggve_h2_analysis(Y, X):
	X_no_high_ld_snps = remove_high_ld_snps(X, .999999999)
	#g_var, e_var = full_inference(np.transpose(X_no_high_ld_snps), Y)
	g_var = full_inference_fixed_residual_var(np.transpose(X_no_high_ld_snps), Y, 1.0)
	return g_var



# Command line args
gtex_genotype_dir = sys.argv[1]
gtex_expression_dir = sys.argv[2]
gtex_covariate_dir = sys.argv[3]
cis_window = sys.argv[4]
tissue_name = sys.argv[5]
chrom_num = sys.argv[6]
gene_h2_estimate_dir = sys.argv[7]


# create output root stem
output_root = gene_h2_estimate_dir + 'h2_estimates_' + tissue_name + '_chr' + chrom_num + '_window_' + cis_window 

# Load in expression data
expr_mat, expr_ind_ids, expr_gene_ids, expr_gene_tss = load_in_expression_data(gtex_expression_dir, tissue_name, chrom_num)

# Load in Covariate data
cov_mat, cov_names, cov_ind_ids = load_in_covariate_data(gtex_covariate_dir, tissue_name)

# Load in genotype data
genotype_stem = gtex_genotype_dir + tissue_name + '_GTEx_v8_genotype_EUR_' + chrom_num
genotype_data = load_in_genotype_data(gtex_genotype_dir, tissue_name, chrom_num)

# Quick error checks to make sure individuals line up between expression, covariates, and genotype
if np.array_equal(cov_ind_ids, expr_ind_ids) == False or np.array_equal(expr_ind_ids, genotype_data['sample_id']) == False:
	print('assumption eroror')
	pdb.set_trace()


##################################
# Run heritability analysis
##################################
# loop through genes and run heritability analysis seperately in each gene
n_genes = len(expr_gene_ids)
aa = []
bb = []
for gene_iter in range(n_genes):
	# Extract relevent information for this gene
	gene_id = expr_gene_ids[gene_iter]
	gene_tss = expr_gene_tss[gene_iter]
	expr_vec = expr_mat[gene_iter, :]

	# Residualize gene expression
	resid_expr_vec = residualize_gene_expression(expr_vec, cov_mat)

	# Extract array of snps in cis window of gene
	cis_snps = (genotype_data['position'] < gene_tss + int(cis_window)) & (genotype_data['position'] >= gene_tss - int(cis_window))
	n_cis_snps = np.sum(cis_snps)
	if n_cis_snps < 20:
		continue

	# RUN GGVE variance estimation
	raw_cis_genotype = genotype_data['G'][:, cis_snps]
	standardized_cis_genotype = standardize_columns_of_matrix(raw_cis_genotype)
	ggve_h2 = run_ggve_h2_analysis(resid_expr_vec, standardized_cis_genotype)


	# RUN GREML ANALYSIS
	tmp_output_stem = output_root + '_tmp_results_' + gene_id + '_'
	#greml_h2, greml_h2_se, greml_h2_p = run_greml_h2_analysis(chrom_num, gene_tss, expr_vec, genotype_stem, cov_mat, cis_window, expr_ind_ids, tmp_output_stem)

	# RUN GREML ANALYSIS WITH NO COVARIATES
	tmp_output_stem = output_root + '_tmp_results_nc_' + gene_id + '_'
	greml_nc_h2, greml_nc_h2_se, greml_nc_h2_p = run_greml_no_covariate_h2_analysis(chrom_num, gene_tss, resid_expr_vec, genotype_stem, cis_window, expr_ind_ids, tmp_output_stem)


	aa.append(ggve_h2)
	bb.append(greml_nc_h2)
	print(aa)
	print(bb)
	print(np.corrcoef(aa,bb))
	if gene_iter > 50:
		break
pdb.set_trace()








