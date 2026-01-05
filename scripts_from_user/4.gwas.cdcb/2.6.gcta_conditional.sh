cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip

gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64


$gcta --bfile ../9.emmax/cohort.chr_1.snp.maf5.thin10000 --make-grm --out cohort.chr_1.snp.maf5.thin10000.grm


n_sample=173
awk -v n_sample=173 'FNR==NR{a[$2]=$2" "$4" "$3" "$6;next}{if(a[$1]){print a[$1],$2,$3,$4,n_sample}}' \
    ../9.emmax/All2_BTAN_WGS.var.bfile.maf5.allele_freq.afreq \
    ../9.emmax/GWAS173_PG.fat_var.emx_k.ps > GWAS173_PG.fat_var.cojo.in

cat ../9.emmax/sig.fat_sv.emx_k.sig.ps | cut -f 1 | while read id ;
do
    let "n_run+=1"
    echo $id > cond.snplist.$n_run
    $gcta --bfile ../9.emmax/All2_BTAN_WGS.var.bfile.maf5 --maf 0.05 --cojo-file GWAS173_PG.fat_var.cojo.in --cojo-cond cond.snplist.$n_run --cojo-p 3.950393e-09 --out cojo.var.$n_run &
done

# $gcta --bfile ../9.emmax/All2_BTAN_WGS.var.bfile.maf5 --chr 14 --maf 0.05 --cojo-file GWAS173_PG.fat_var.cojo.in --cojo-cond cond.snplist --cojo-p 3.950393e-09 --out test_chr14


# independently associated SNPs
$gcta --bfile ../9.emmax/All2_BTAN_WGS.var.bfile.maf5 --maf 0.05 --cojo-file GWAS173_PG.fat_var.cojo.in --cojo-slct --cojo-p 3.950393e-09 --out cojo.slct

$gcta --bfile ../9.emmax/All2_BTAN_WGS.var.bfile.maf5 --maf 0.05 --cojo-file GWAS173_PG.fat_var.cojo.in --cojo-top-SNPs 100 --out cojo.top100

cat test_chr14.top100.cma.cojo | grep -v "SNP" | grep -v "s" | wc -l


#  a fixed number of independently associated SNPs without a p-value threshold
--cojo-top-SNPs 10\


