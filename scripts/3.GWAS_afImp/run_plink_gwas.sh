# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${PHENO_DIR:?PHENO_DIR is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
: "${SLURM_PARTITION:?SLURM_PARTITION is unset - see config.sh.example at the repository root}"
# --------------------------

#cat > ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6

chrom=29
nthreads=$SLURM_CPUS_PER_TASK

if [[ $cov == "F" ]]; then
    plink2 --threads $nthreads --bfile $pre_gdata \
        --chr-set $chrom \
        --pheno iid-only $ph --pheno-name $trait \
        --no-categorical --glm allow-no-covars \
        --no-input-missing-phenotype \
        --out $out_pre.pc$pc
else
    plink2 --threads $nthreads --bfile $pre_gdata --chr-set 29 \
        --pheno iid-only $ph --pheno-name $trait \
        --covar $cov --covar-name $(printf "PC%d " $(seq 1 $pc)) \
        --no-categorical --glm hide-covar \
        --no-input-missing-phenotype \
        --out $out_pre.pc$pc
fi
flt=$(awk '$0 ~ "variants loaded from"{print $1}' $out_pre.pc$pc.log)
awk '$15 < 0.05/flt && $10 == "ADD"{print $3, $7, $4, $9, $12, $13, $15, $11}' \
    flt=$flt $out_pre.pc$pc.$trait.glm.linear > $out_pre.pc$pc.$trait.summary_stats.txt
cat $out_pre.pc$pc.$trait.summary_stats.txt | 
    awk '$1 !~ "SNP" && $1 !~ "s"{print $1}' > $out_pre.pc$pc.$trait.svlist
cat $out_pre.pc$pc.$trait.summary_stats.txt | 
    awk '{split($1,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc.$trait /" \
    > $out_pre.pc$pc.$trait.svlist.count
EOF


nthreads=16
for var in sv var gsnp; do

pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
out_pre=cdcb_gwas_nocv.$var
trait=Milk
cov=All2_BTAN_WGS.var.bfile.maf1.eigenvec
pc="2"
ph=ph.hol.traits

for trait in fat milk pro; do

sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done
done

nthreads=16
for var in sv var gsnp; do

pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
out_pre=cdcb_gwas_nocv.$var
trait=Milk
cov=F
pc=F
ph=cdcb.hol.traits

for trait in fat milk pro; do

sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done
done



cd ${PROJECT_ROOT}/hd_gwas/genotypes/1.svImp_wgs

nthreads=48
pre_gdata=umd50ksWGS_HolrePan.filter
out_pre=cdcb-umd50k_gwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=2
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done


cd ${PROJECT_ROOT}/hd_gwas/genotypes/4.snp.gwas
nthreads=48
pre_gdata=../allseq_1kbulls.hol.chrall.maf01
out_pre=cdcb-umd50k_snpgwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=10
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done


cd ${PROJECT_ROOT}/hd_gwas/genotypes/2.svImp_ld

nthreads=48
pre_gdata=umd50ksLD_HolrePan.filter2
out_pre=cdcb-umd50k_ldgwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=2
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done



cd ${PROJECT_ROOT}/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=10
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

done


cd ${PROJECT_ROOT}/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=10
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A ${SLURM_ACCOUNT} \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

done


cd ${PROJECT_ROOT}/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=${PHENO_DIR}/hol.chip.maf01.eigenvec
pc=10
ph=${PHENO_DIR}/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do

sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done

