
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64
sbatch -A bull_scr \
        --job-name=cdcb-umd50k_gcta \
        --cpus-per-task=48 \
        --error=cdcb-umd50k_gcta.err \
        --output=cdcb-umd50k_gcta.out \
        --time "1-00:00:00" \
        --wrap="
$gcta --bfile hol.chip --make-grm --autosome-num 29 \
    --thread-num 48 --out hol.chip.grm
$gcta --grm hol.chip.grm --make-bK-sparse 0.05 --out hol.chip.sp_grm
"

cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_grm.sh <<'EOF'
pre_gdata=$1
nall=$2
nc=$3

chrom=29
nthreads=$SLURM_CPUS_PER_TASK
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64
mkdir -p $pre_gdata.grm
$gcta --bfile $pre_gdata --make-grm-part $nall $nc --autosome-num $chrom --chr 1 --maf 0.1  --max-maf 0.4 \
    --thread-num $nthreads --out $pre_gdata.grm/$pre_gdata.grm
EOF

pre_gdata=umd50ksWGS_HolrePan.filter
nthreads=16
nall=100
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs
mkdir -p $pre_gdata.grm

seq 1 $nall | while read ii; do
    sbatch -A bull_scr \
        --job-name=$ii.grm \
        --cpus-per-task=$nthreads \
        --error=$pre_gdata.grm/$ii.grm.err \
        --output=$pre_gdata.grm/$ii.grm.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_grm.sh \
    $pre_gdata $nall $ii
"
done

cat $pre_gdata.grm/$pre_gdata.grm.part_${nall}_*.grm.id > $pre_gdata.grm.id
cat $pre_gdata.grm/$pre_gdata.grm.part_${nall}_*.grm.bin > $pre_gdata.grm.bin
cat $pre_gdata.grm/$pre_gdata.grm.part_${nall}_*.grm.N.bin > $pre_gdata.grm.N.bin


#cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6
grm=$7

chrom=29
nthreads=$SLURM_CPUS_PER_TASK
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

if [[ $cov == "F" ]]; then
    $gcta --fastGWA-mlm --bfile $pre_gdata --grm-sparse $grm \
        --pheno $ph --autosome-num $chrom \
        --out $out_pre --thread-num $nthreads
else
    $gcta --fastGWA-mlm --bfile $pre_gdata --grm-sparse $grm \
        --pheno $ph --autosome-num $chrom --qcovar $cov \
        --out $out_pre.pc$pc --thread-num $nthreads
fi
EOF


#cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6
grm=$7

chrom=29
nthreads=$SLURM_CPUS_PER_TASK
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

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

#### 1.svImp_wgs ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs

nthreads=16
pre_gdata=umd50ksWGS_HolrePan.filter
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.sp_grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait
out_pre=cdcb-umd50k_gctaPc10fastmlm_$trait

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

out_pre=cdcb-umd50k_gctafastmlm_$trait
cov=F
pc=F
sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done

#### 2.svImp_ld ==================================================================================================================================
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld

nthreads=16
pre_gdata=umd50ksLD_HolrePan.filter2
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.sp_grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait
out_pre=cdcb-umd50k_gctaPc10fastmlm_$trait

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

out_pre=cdcb-umd50k_gctafastmlm_$trait
cov=F
pc=F
sbatch -A bull_scr --partition bigmem \
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

#### 5.svImp_svwgs ==================================================================================================================================


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs

nthreads=16
pre_gdata=umd50ksWGSsv_HolrePan.filter
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.sp_grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait
out_pre=cdcb-umd50ksv_gctaPc10fastmlm_$trait

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

out_pre=cdcb-umd50ksv_gctafastmlm_$trait
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


#### 6.svImp_svld ==================================================================================================================================


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/6.svImp_svld

nthreads=16
pre_gdata=umd50ksLDsv_HolrePan.filter
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.sp_grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait
out_pre=cdcb-umd50kldsv_gctaPc10fastmlm_$trait

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "4-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

out_pre=cdcb-umd50kldsv_gctafastmlm_$trait
cov=F
pc=F
sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "4-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_fastmlm.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done
done





