# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${PHENO_DIR:?PHENO_DIR is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${PHENO_DIR}

 awk  'BEGIN{FS=OFS=","} FNR==NR{a[$1]=$2","$3","$4","$5;next} {print $0,a[$1]}' \
        ${PROJECT_ROOT}/hd_gwas/hol.type.pheno.csv \
        ${PROJECT_ROOT}/hd_gwas/hol.yld.pheno.csv > \
        ${PROJECT_ROOT}/hd_gwas/hol.yld_type.pheno.csv

cat ${PROJECT_ROOT}/hd_gwas/hol.yld.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        ${PHENO_DIR}/hol.yld.pheno.phe2

cat ${PROJECT_ROOT}/hd_gwas/hol.type.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        ${PHENO_DIR}/hol.type.pheno.phe2

cat ${PROJECT_ROOT}/hd_gwas/hol.yld_type.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        ${PHENO_DIR}/hol.yld_type.pheno.phe2

cat ${PHENO_DIR}/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$2}' \
    > ${PHENO_DIR}/hol.yld.pheno.Milk
cat ${PHENO_DIR}/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$4}' \
    > ${PHENO_DIR}/hol.yld.pheno.Fat
cat ${PHENO_DIR}/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$6}' \
    > ${PHENO_DIR}/hol.yld.pheno.Protein

cat ${PHENO_DIR}/hol.type.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$2}' | awk '$3 != ""' \
    > ${PHENO_DIR}/hol.yld.pheno.Stature

cat ${PHENO_DIR}/hol.type.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$4}' | awk '$3 != ""' \
    > ${PHENO_DIR}/hol.yld.pheno.Body_depth


cat ${PHENO_DIR}/hol.yld.pheno.phe2 | 
    sed -e '1d' -e 's/^/0\t/g' > ${PHENO_DIR}/gcta.hol.phen
cat ${PHENO_DIR}/hol.chip.maf01.eigenvec | 
    sed -e '1d' -e 's/^/0\t/g' > ${PHENO_DIR}/gcta.hol.10PCs


plink2 --threads 8 --vcf ../hol.chip.maf01.vcf.gz --chr-set 29 --const-fid --pca 10 --out hol.chip.maf01

#cd ${PROJECT_ROOT}/hd_gwas/genotypes/2.svImp_ld
cd ${PROJECT_ROOT}/hd_gwas/genotypes/1.svImp_wgs

trait=Milk
#pre_gdata=umd50ksLD_HolrePan.filter2
pre_gdata=umd50ksWGS_HolrePan.filter
out_pre=cdcb-umd50k_gwas_pc10
pc=10
nthreads=48
phe=${PHENO_DIR}/

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=$trait.GWAS \
        --cpus-per-task=$nthreads \
        --error=$trait.GWAS.err \
        --output=$trait.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
plink2 --threads $nthreads --bfile $pre_gdata --chr-set 29 \
       --pheno iid-only $phe/hol.yld.pheno.phe2 --pheno-name $trait \
       --covar $phe/hol.chip.maf01.eigenvec --covar-name $(printf "PC%d " $(seq 1 $pc)) \
        --no-categorical --glm --out $out_pre
flt=\$(awk '\$0 ~ \"variants loaded from\"{print \$1}' $out_pre.log)
awk '\$15 < 0.05/flt && \$10 == \"ADD\"{print \$3, \$7, \$4, \$9, \$12, \$13, \$15, \$11}' \
       flt=\$flt $out_pre.$trait.glm.linear > $out_pre.$trait.summary_stats.txt
cat $out_pre.$trait.summary_stats.txt |
        awk '\$1 !~ \"SNP\" && \$1 !~ \"s\"{print \$1}' > $out_pre.$trait.svlist
"
done

