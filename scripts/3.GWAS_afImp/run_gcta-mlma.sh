# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${PHENO_DIR:?PHENO_DIR is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
: "${SLURM_PARTITION:?SLURM_PARTITION is unset - see config.sh.example at the repository root}"
# --------------------------

#cat > ${PROJECT_ROOT}/hd_gwas/genotypes/gcta_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6
grm=$7

chrom=29
nthreads=$SLURM_CPUS_PER_TASK
gcta=${PROJECT_ROOT}/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

if [[ $cov == "F" ]]; then
    $gcta --mlma --bfile $pre_gdata --grm $grm \
        --pheno $ph --autosome-num $chrom \
        --out $out_pre.pc$pc --thread-num $nthreads
else
    $gcta --mlma --bfile $pre_gdata --grm $grm \
        --pheno $ph --autosome-num $chrom --qcovar $cov \
        --out $out_pre.pc$pc --thread-num $nthreads
fi
EOF

#### 5.svImp_svwgs ==================================================================================================================================


cd ${PROJECT_ROOT}/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=${PHENO_DIR}/gcta.hol.10PCs
pc=10
grm=${PROJECT_ROOT}/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=${PHENO_DIR}/hol.yld.pheno.$trait
out_pre=cdcb-umd50ksv_gctamlma_$trait
cov=F
pc=F

sbatch -A ${SLURM_ACCOUNT} --partition bigmem \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=30G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "14-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

done


#### 6.svImp_svld ==================================================================================================================================

cd ${PROJECT_ROOT}/hd_gwas/genotypes/6.svImp_svld

nthreads=48
pre_gdata=umd50ksLDsv_HolrePan.filter
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=${PHENO_DIR}/gcta.hol.10PCs
pc=10
grm=${PROJECT_ROOT}/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=${PHENO_DIR}/hol.yld.pheno.$trait
out_pre=cdcb-umd50kldsv_gctamlma_$trait
cov=F
pc=F
sbatch -A ${SLURM_ACCOUNT} --partition bigmem \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=30G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "14-00:00:00" \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done


#### 1.svImp_wgs ==================================================================================================================================

cd ${PROJECT_ROOT}/hd_gwas/genotypes/1.svImp_wgs

nthreads=16
pre_gdata=umd50ksWGS_HolrePan.filter
grm=${PROJECT_ROOT}/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=${PHENO_DIR}/hol.yld.pheno.$trait
out_pre=cdcb-umd50k_gctamlma_$trait
cov=F
pc=F
sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done

#### 2.svImp_ld ==================================================================================================================================
cd ${PROJECT_ROOT}/hd_gwas/genotypes/2.svImp_ld

nthreads=16
pre_gdata=umd50ksLD_HolrePan.filter2
grm=${PROJECT_ROOT}/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=${PHENO_DIR}/hol.yld.pheno.$trait
out_pre=cdcb-umd50k_gctamlma_$trait
cov=F
pc=F
sbatch -A ${SLURM_ACCOUNT} --partition ${SLURM_PARTITION} \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --wrap="
bash ${PROJECT_ROOT}/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done


