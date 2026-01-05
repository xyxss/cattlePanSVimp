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

#### 5.svImp_svwgs ==================================================================================================================================


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/00.svImp_r85sv

plink2 --threads 8 --bfile ../umdsample.svwgssnp.afterImput.filter \
    --extract ../umd50ksWGS_HolrePan.filter.r85.sv --chr-set 29 --make-pgen --out umd50ksWGSsv_HolrePan.filter.r85


plink2 --threads 8 --bfile ../umdsample.svwgssnp.afterImput.filter \
    --extract ../umd50ksWGS_HolrePan.filter.r85.sv --chr-set 29 --make-bed --out umd50ksWGSsv_HolrePan.filter.r85
    

nthreads=16
pre_gdata=umd50ksWGSsv_HolrePan.filter.r85
#out_pre=cdcb-umd50k_gcta2_Fat
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein Stature Body_depth; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.$trait

out_pre=cdcb-umd50ksv_gctamlma_$trait

cov=F
pc=F

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}.GCTA \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=20G \
        --error=${out_pre}.GCTA.err \
        --output=${out_pre}.GCTA.out \
        --time "14-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/gcta_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

done


65151



nthreads=16
pre_gdata=umd50ksWGSsv_HolrePan.filter.r85


#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cat $pre_gdata.pvar | awk 'BEGIN{print "SNP"} $0 !~ "#"{print $3}' > $pre_gdata.snp_info.csv

grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm
cov=F
pc=F
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}-pc${pc}.slemm \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=8G \
        --error=logs/${out_pre}-pc${pc}.slemm.err \
        --output=logs/${out_pre}-pc${pc}.slemm.out \
        --time "14-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/slemm_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

done

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait

cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt

awk 'NR == 1 || $9 < 0.05/65151' $out_pre.pc$pc.gwa.txt \
    >  $out_pre.pc$pc.gwa.sig
cat $out_pre.pc$pc.gwa.sig | 
    awk '$3 !~ "SNP" && $3 !~ "s"{print}' > $out_pre.pc$pc.svlist
cat $out_pre.pc$pc.gwa.sig | 
    awk '{split($3,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc /" \
    > $out_pre.pc$pc.svlist.count
done

