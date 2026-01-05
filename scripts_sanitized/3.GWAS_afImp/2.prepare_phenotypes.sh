


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes

 awk  'BEGIN{FS=OFS=","} FNR==NR{a[$1]=$2","$3","$4","$5;next} {print $0,a[$1]}' \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.type.pheno.csv \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld.pheno.csv > \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.type.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.type.pheno.phe2

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv | 
        sed -e 's/,/\t/g' -e 's/1_HOL_IID/#IID/' > \
        /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld_type.pheno.phe2

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$2}' \
    > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.Milk
cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$4}' \
    > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.Fat
cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$6}' \
    > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.Protein

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.type.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$2}' | awk '$3 != ""' \
    > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.Stature

cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.type.pheno.phe2 | 
    awk 'BEGIN{FS=OFS="\t"} NR > 1{print 0"\t"$1"\t"$4}' | awk '$3 != ""' \
    > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.Body_depth


cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2 | 
    sed -e '1d' -e 's/^/0\t/g' > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.phen
cat /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec | 
    sed -e '1d' -e 's/^/0\t/g' > /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/gcta.hol.10PCs


plink2 --threads 8 --vcf ../hol.chip.maf01.vcf.gz --chr-set 29 --const-fid --pca 10 --out hol.chip.maf01

#cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs

trait=Milk
#pre_gdata=umd50ksLD_HolrePan.filter2
pre_gdata=umd50ksWGS_HolrePan.filter
out_pre=cdcb-umd50k_gwas_pc10
pc=10
nthreads=48
phe=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
sbatch -A bull_scr \
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

