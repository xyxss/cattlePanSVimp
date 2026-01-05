


/90daydata/bull_age/liu.yang/imputation/hd_gwas



#ref https://groups.google.com/a/soe.ucsc.edu/g/genome/c/pmpGHtsyiP4


Here are the details for the two assemblies:
GCF_002263795.3_ARS-UCD2.0 cattle (Hereford L1 Dominette 01449 42190680 v2.0 2023 USDA)
GCF_000003055.6_Bos_taurus_UMD_3.1.1 cattle (Hereford UMD_3.1.1 2014)

The LiftOver files are available here:
https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/liftOver/GCF_002263795.3ToGCF_000003055.6.over.chain.gz
https://hgdownload.soe.ucsc.edu/hubs/GCF/000/003/055/GCF_000003055.6/liftOver/GCF_000003055.6ToGCF_002263795.3.over.chain.gz

And you can view these genomes and the chain tracks (that display these alignments) here:
https://genome.ucsc.edu/h/GCF_002263795.3
https://genome.ucsc.edu/h/GCF_000003055.6

I hope this is helpful. If you have any further questions, please reply to gen...@soe.ucsc.edu. All messages sent to that address are archived on a publicly-accessible Google Groups forum. If your question includes sensitive data, you may send it instead to genom...@soe.ucsc.edu.

wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/003/055/GCF_000003055.6/liftOver/GCF_000003055.6ToGCF_002263795.3.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/liftOver/GCF_002263795.3ToBosTau9.over.chain.gz


awk 'FNR==NR{a[$1]=$2;next} {split($2,aa,":");print a[aa[1]]"\t"aa[2]"\t"aa[2]}' /90daydata/bull_age/liu.yang/ref/GCF_002263795.3_ARS-UCD2.0.chr.tr hol.chip.bim > hol.chip.bed
#775032

liftOver hol.chip.bed GCF_000003055.6ToGCF_002263795.3.over.chain ARS-UCD2.0.bed unMapped.bed

#ARS-UCD2.0.bed
#774316

awk 'FNR==NR{a[$2]=$1;next} {print a[$1]"\t"$2"\t"$2}' /90daydata/bull_age/liu.yang/ref/GCF_002263795.3_ARS-UCD2.0.chr.tr ARS-UCD2.0.bed  > GCF_002263795.3_ARS-UCD2.0.bed

cat GCF_002263795.3_ARS-UCD2.0.bed | awk '{print $1"\t"$2}' > GCF_002263795.3_ARS-UCD2.0.pos





plink2 --pfile filtered_data --pheno phenotype.txt --covar covariates.txt --glm --out gwas_results

--pheno phenotype.txt
FID IID PHENOTYPE
1    1   1
1    2   2

--covar covariates.txt

FID IID AGE SEX PC1 PC2 PC3
1    1   45  1   0.1 0.2 0.3
1    2   32  2   0.2 0.3 0.4



awk '$12 < 5e-8' gwas_results.PHENO1.glm.logistic > significant_results.txt

# Convert to PLINK 2 format

# Perform QC
plink2 --pfile step1_data --mind 0.1 --geno 0.05 --maf 0.01 --make-pgen --out qc_data

# Calculate PCs
plink2 --pfile qc_data --pca 10 --out pca_results

# Run GWAS with covariates
plink2 --pfile qc_data --pheno phenotype.txt --covar covariates.txt --glm --out gwas_results



#### phepotype to plink
awk 'FNR==NR{a[$1]++;next} {split($1,aa,",");if(a[aa[1]])print} ' allseq_1kbulls.hol.chrall.psam ../hol.yld.pheno.csv  | wc -l
50309

awk 'FNR==NR{a[$1]++;next} {split($1,aa,",");if(a[aa[1]])print} ' allseq_1kbulls.hol.chrall.psam ../hol.yld.pheno.csv

cat ../hol.yld.pheno.csv | awk -f phcor.awk > hol.yld.pheno.phe


plink2 --threads 8 --pfile allseq_1kbulls.hol.chrall --chr-set 29 --const-fid --pca 10 --out allseq_1kbulls.hol.chrall.pca


plink2 --pfile qc_data --pheno iid-only phenotype.txt --pheno-name $trait --covar covariates.txt --glm --out gwas_results

ls /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/*gz | 
    while read id; do
        sbatch -A bull_scr -c 2 --wrap="
bcftools stats --threads 2 $id > ${id/vcf.gz/vcfstats}
"
    done





####
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
...


plink2 --pfile input_prefix --make-bed --out gcta_input

#$gcta --bfile ../9.emmax/cohort.chr_1.snp.maf5.thin10000 --make-grm --out cohort.chr_1.snp.maf5.thin10000.grm

# Conditional Analysis with GCTA
gcta64 --bfile gcta_input \
       --cojo-file summary_stats.txt \
       --cojo-cond SNP1 SNP2 ... \
       --out gcta_conditional

# Generate Summary Statistics File
plink2 --pfile input_prefix \
       --glm cols=+a1freq,+err \
       --out glm_results

# Run Joint Analysis (Optional)
gcta64 --bfile gcta_input \
       --cojo-file summary_stats.txt \
       --cojo-joint \
       --out gcta_joint

Explanation of Columns
Column Name	Description
#CHROM	Chromosome number of the SNP.
POS	Physical position of the SNP on the chromosome (in base pairs).
ID	SNP identifier (e.g., rsID or ., if no identifier is assigned).
REF	Reference allele at this position in the genome.
ALT	Alternate allele (the tested allele in GWAS).
PROVISIONAL_REF?	Indicates whether the reference allele was "provisionally" set during the analysis.
A1	Effect allele (same as ALT if not otherwise specified).
OMITTED	Placeholder for dropped SNPs or issues during computation (e.g., no variation).
A1_FREQ	Frequency of the effect allele (A1) in the sample.
TEST	The type of association test performed (e.g., ADD for additive, DOM for dominant, etc.).
OBS_CT	Number of observed samples for this SNP.
BETA	Effect size (linear regression coefficient).
SE	Standard error of the effect size.
T_STAT	T-statistic for the test.
P	P-value for the test of association.
ERRCODE	Error code (non-zero if there were issues for this SNP).

