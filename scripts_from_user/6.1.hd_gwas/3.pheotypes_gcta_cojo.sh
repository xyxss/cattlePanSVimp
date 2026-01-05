


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes


### conditional analysis

cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip

plink2 --threads 8 --pfile ../genotypes/umdsample.svwgssnp.afterImput.filter --make-bed --out ../genotypes/umdsample.svwgssnp.afterImput.filter
#$gcta --bfile ../9.emmax/cohort.chr_1.snp.maf5.thin10000 --make-grm --out cohort.chr_1.snp.maf5.thin10000.grm


trait=Milk
nthreads=48
#pre_gdata=umd50ksLD_HolrePan.filter2
pre_gdata=umd50ksWGS_HolrePan.filter

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
flt=$(awk '$0 ~ "variants loaded from"{print $1}' $pre_gdata.gwas.log)
awk '$15 < 0.05/flt && $10 == "ADD"{print $3, $7, $4, $9, $12, $13, $15, $11}' \
       flt=$flt $pre_gdata.gwas.$trait.glm.linear > $pre_gdata.gwas.$trait.summary_stats.txt
cat $pre_gdata.gwas.$trait.summary_stats.txt | 
       awk '$1 !~ "SNP" && $1 !~ "s"{print $1}' > $pre_gdata.gwas.$trait.svlist
gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64
cat > gcta_conditional.$trait.sh <<EOF
$gcta  --threads $nthreads --bfile $pre_gdata \
       --cojo-file $pre_gdata.gwas.$trait.summary_stats.txt \
       --cojo-slct \
       --out cojo-$pre_gdata.$trait.gcta_slct
$gcta  --threads $nthreads --bfile $pre_gdata \
       --cojo-file $pre_gdata.gwas.$trait.summary_stats.txt \
       --cojo-joint \
       --out cojo-$pre_gdata.$trait.gcta_joint
$gcta  --threads $nthreads --bfile $pre_gdata \
       --cojo-file $pre_gdata.gwas.$trait.summary_stats.txt \
       --cojo-top-SNPs 100 \
       --out cojo-$pre_gdata.$trait.gcta_top100
$gcta  --threads $nthreads --bfile $pre_gdata \
       --cojo-file $pre_gdata.gwas.$trait.summary_stats.txt \
       --extract $pre_gdata.$trait.svlist \
       --cojo-joint \
       --out cojo-$pre_gdata.$trait.gcta_top10
EOF
cat gcta_conditional.$trait.sh | 
       while read line; do
       sbatch -A bull_scr \
        --job-name=$trait.gcta_conditional \
        --cpus-per-task=$nthreads \
        --error=gcta_conditional.$trait.err \
        --output=gcta_conditional.$trait.out \
        --time "1-00:00:00" \
        --wrap="
$line
"
done 
done

