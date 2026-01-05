cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85

#### cojo freq
awk -F, 'NR==1{print "#IID\tMilk"; next} {print $1"\t"$2}' /90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv \
    > ph.hol.milk

trait=Milk
nthreads=48
#pre_gdata=umd50ksLD_HolrePan.filter2
pre_gdata=umd50ksWGS_HolrePan.filter.r85

out_pre=cdcb_gwas_nocv.milk
trait=Milk
ph=ph.hol.milk

plink2 --threads $nthreads --bfile $pre_gdata \
    --chr-set 29 \
    --pheno iid-only $ph --pheno-name $trait \
    --no-categorical --glm allow-no-covars \
    --no-input-missing-phenotype \
    --out $out_pre

####

awk 'BEGIN{print "SNP A1 A2 Freq b se p"} 
    FNR==NR{if($10 == "ADD" && $1 != "ID")a[$3]=$9;next} 
    {print $3,$4,$5,a[$3],$6,$7,$9,"50299"}' \
    cdcb_gwas_nocv.milk.Milk.glm.linear \
    <(sed '1d' cdcb-umd50ksv_slemm_Stature.pcF.gwa.chr3.txt) \
    > cdcb-umd50ksv_slemm_Stature.pcF.gwa.chr3.cojo.in

trait=Milk
nthreads=48
#pre_gdata=umd50ksLD_HolrePan.filter2
pre_gdata=umd50ksWGS_HolrePan.filter.r85
flt=$(awk '$0 ~ "variants loaded from"{print $1}' $pre_gdata.log)

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
0.05/6048794
awk '$9 < 0.05/flt{print $3, $7, $4, $9, $12, $13, $15, $11}' \
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



gcta=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/10.gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

flt=$(awk '$0 ~ "variants loaded from"{print $1}' $pre_gdata.log)

awk '$15 < 0.05/flt && $10 == "ADD"{print $3, $7, $4, $9, $12, $13, $15, $11}' \
       flt=$flt cdcb_gwas_nocv.milk.Milk.glm.linear > cdcb_gwas_nocv.milk.Milk.summary_stats.txt


trait=Stature
chr=3
out_pre=cdcb-umd50ksv_slemm_$trait

pre_gdata=umd50ksWGS_HolrePan.filter.r85
$gcta  --threads $nthreads --bfile $pre_gdata --chr $chr \
       --cojo-file cdcb_gwas_nocv.milk.Milk.summary_stats.txt \
       --cojo-slct \
       --out test2.cojo
