cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/9.emmax


plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-snps.vcf.gz --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.05 --make-bed --out All2_BTAN_WGS.var-snp.bfile.maf5


shuf -n 61249 All2_BTAN_WGS.var-snp.bfile.maf.bim | cut -f 2 > All2_BTAN_WGS.var-snp.bfile.maf.snplist

plink2 --threads $nthreads --bfile All2_BTAN_WGS.var-snp.bfile.maf --extract All2_BTAN_WGS.var-snp.bfile.maf.snplist --const-fid 0 --chr-set 29 --chr 1-29  --make-bed --out All2_BTAN_WGS.var-snp-dsv.bfile.maf


for typ in var-snp-dsv var-snp; do
plink --threads $nthreads --bfile All2_BTAN_WGS.$typ.bfile.maf --recode12 --output-missing-genotype 0 --transpose --out All2_BTAN_WGS.$typ.bfile.maf5.t --chr-set 29 
done



trait="fat"
nt=8
for typ in var-snp-dsv var-snp; do
    ls ind_154animals_dPTA.*.txt | cut -d. -f2 | while read trait; do
    echo "./emmax-intel64 -v -d 10 -t All2_BTAN_WGS.$typ.bfile.maf5.t -p ind_154animals_dPTA.$trait.txt -k cohort.chr_1.snp.maf5.thin10000.t.aIBS.kinf -o GWAS173_PG.${trait}_$typ.emx_k"
done 
done > run_line1.txt

cat run_line1.txt | while read line; do
let "n_run+=1"
    sbatch -A bull_scr \
            -D $PWD \
        --export=ALL \
        -J gwas_$n_run \
        -c $nt \
        -o logs/gwas_$n_run.out \
        --wrap="
$line
        "
done
