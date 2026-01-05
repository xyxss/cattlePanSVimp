cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/4.snp.gwas

nthreads=48
pre_gdata=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/allseq_1kbulls.hol.chrall.maf01
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait
out_pre=cdcb-umd50k_gcta2_$trait

sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

out_pre=cdcb-umd50k_gcta_$trait
cov=F
pc=F
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done
done
