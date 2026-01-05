cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/plink_gwas


for var in sv var gsnp; do
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64
sbatch -A bull_scr \
        --job-name=cdcb-umd50k_gcta \
        --cpus-per-task=48 \
        --error=cdcb-umd50k_gcta.err \
        --output=cdcb-umd50k_gcta.out \
        --time "2-00:00:00" \
        --wrap="
$gcta --bfile All2_BTAN_WGS.$var.bfile.maf1 --make-grm --autosome-num 29 \
    --thread-num 48 --out All2_BTAN_WGS.$var.grm
$gcta --grm All2_BTAN_WGS.$var.grm --make-bK-sparse 0.05 --out All2_BTAN_WGS.$var.sp_grm
"
done

for var in sv var gsnp; do
nthreads=16
pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
grm=All2_BTAN_WGS.gsnp.sp_grm

for trait in milk fat pro; do
ph=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/9.emmax/ind_154animals_dPTA.$trait.txt

out_pre=cdcb-173wg-${var}_gctafastmlm_$trait
cov=F
pc=F
sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "5-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done
done


for var in sv var gsnp; do
for trait in milk fat pro; do
out_pre=cdcb-173wg-${var}_gctafastmlm_$trait
flt=$(awk '$1 ~ "Saved"{print $2}' $out_pre.log)
awk '$10 < 0.05/flt' \
    flt=$flt $out_pre.fastGWA > $out_pre.summary_stats.txt
cat $out_pre.summary_stats.txt | 
    awk '$2 !~ "SNP" && $2 !~ "s"{print $1}' > $out_pre.svlist
cat $out_pre.summary_stats.txt | 
    awk '{split($2,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre /" \
    > $out_pre.svlist.count
done
done

