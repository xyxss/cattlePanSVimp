
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/sv_gwas

plink2 --threads 8 --pfile allseq_1kbulls.hol.chrall.maf01 --chr-set 29 --const-fid --pca 10 --out allseq_1kbulls.hol.chrall.pca

plink2 --pfile qc_data --pheno iid-only phenotype.txt --pheno-name $trait --covar covariates.txt --glm --out gwas_results