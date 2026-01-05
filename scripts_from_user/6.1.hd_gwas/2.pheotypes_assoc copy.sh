


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

awk '{print $1,$2,$3}' hol.chip.maf01.eigenvec


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



cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs

 wc -l *.summary_stats.txt *svlist
   1984730 cdcb-umd50k_gwas.pc10.Fat_Percent.summary_stats.txt
   4174557 cdcb-umd50k_gwas.pc10.Fat.summary_stats.txt
   5244913 cdcb-umd50k_gwas.pc10.Milk.summary_stats.txt
   2834891 cdcb-umd50k_gwas.pc10.Pro_Percent.summary_stats.txt
   5917420 cdcb-umd50k_gwas.pc10.Protein.summary_stats.txt
      9677 cdcb-umd50k_gwas.pc10.Fat_Percent.svlist
     21586 cdcb-umd50k_gwas.pc10.Fat.svlist
     26756 cdcb-umd50k_gwas.pc10.Milk.svlist
     13907 cdcb-umd50k_gwas.pc10.Pro_Percent.svlist
     30413 cdcb-umd50k_gwas.pc10.Protein.svlist

plink2 --threads $nthreads --pfile allseq_1kbulls.hol.chrall.maf01 --chr-set 29 \
        --make-bed --out allseq_1kbulls.hol.chrall.maf01

thr=$(let 4.71807e-09)
cat GWAS173_PG.fat_var-snp.emx_k.ps | 
    awk '
        ($4 < Thr*1e1 && NR % 2 == 1) || 
        ($4 < Thr*1e2 && NR % 5 == 1) || 
        ($4 < Thr*1e3 && NR % 10 == 1) || 
        ($4 > Thr*1e3 && NR % 200 == 1)
        ' \
        Thr=$thr > sig.fat_var-snp.emx_k.base.ps



  14751596
   4174557 cdcb-umd50k_gwas.pc10.Fat.summary_stats.txt
   5244913 cdcb-umd50k_gwas.pc10.Milk.summary_stats.txt
   5917420 cdcb-umd50k_gwas.pc10.Protein.summary_stats.txt
   3657970 cdcb-umd50k_snpgwas.pc10.Fat.summary_stats.txt
   4598731 cdcb-umd50k_snpgwas.pc10.Milk.summary_stats.txt
   5186986 cdcb-umd50k_snpgwas.pc10.Protein.summary_stats.txt
     21586 cdcb-umd50k_gwas.pc10.Fat.svlist
     26756 cdcb-umd50k_gwas.pc10.Milk.svlist
     30413 cdcb-umd50k_gwas.pc10.Protein.svlist

14751596    
cat umd50ksWGS_HolrePan.filter.bim | awk '{split($2,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' 
SNP 11,996,492
DEL 17,984
sDEL 1031465
INS 31,774
sINS 1172459
COMPLEX 29,128
sCOMPLEX 472294