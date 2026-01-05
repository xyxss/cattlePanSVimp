


cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb


## SNP dataset
cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk

bcftools concat --threads $nthreads \
    --allow-overlaps \
    cohort.chr_*.snp.vcf.gz \
    -Oz -o All_BTAN_WGS_gatk.snp.vcf.gz &&
    tabix -f -p vcf All_BTAN_WGS_gatk.snp.vcf.gz

##

## pan SV dataset
code=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/00.codes
cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/

bcftools annotate --threads $nthreads --remove QUAL,FORMAT \
    hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.vcf.gz \
    -Oz -o All2_BTAN_WGS.vcf.gz && 
    tabix -f -p vcf All2_BTAN_WGS.vcf.gz

bcftools view --threads $nthreads -v indels \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz | 
    awk -f  $code/pan.indel.awk |
    bgzip -c > All2_BTAN_WGS-indels.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-indels.vcf.gz &

bcftools view --threads $nthreads -v other,mnp \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz  | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > All2_BTAN_WGS-complex.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-complex.vcf.gz &
wait

bcftools concat --threads $nthreads \
    --allow-overlaps \
    All2_BTAN_WGS-complex.vcf.gz \
    All2_BTAN_WGS-indels.vcf.gz \
    -Oz -o All2_BTAN_WGS-nonsnps.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-nonsnps.vcf.gz

zcat All2_BTAN_WGS-nonsnps.vcf.gz | 
    awk '$1 ~ /#/{print; next}{split($3,a,/:|-|_/);if(a[5] > 50) {print}}' |
    bgzip -c > All2_BTAN_WGS-sv.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-sv.vcf.gz

###

### pan SNP dataset

bcftools view --threads $nthreads -v snps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    All2_BTAN_WGS.vcf.gz  |
    awk -f  $code/pan.snp.awk |
    bgzip -c > All2_BTAN_WGS-snps.vcf.gz && 
    tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz
bcftools sort All2_BTAN_WGS-snps2.vcf.gz -Oz -o All2_BTAN_WGS-snps.vcf.gz
 tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz

bcftools sort All2_BTAN_WGS-snps.vcf.gz -Oz -o All2_BTAN_WGS-snps2.vcf.gz &&
    mv All2_BTAN_WGS-snps2.vcf.gz All2_BTAN_WGS-snps.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-snps.vcf.gz


### pan VAR dataset

bcftools concat --threads 4 \
    --allow-overlaps \
    All2_BTAN_WGS-nonsnps.vcf.gz \
    All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o All2_BTAN_WGS-var.vcf.gz &&
    tabix -f -p vcf All2_BTAN_WGS-var.vcf.gz

### 
### 
# phenotype

cd /90daydata/bull_age/liu.yang/imputation/emmax
cat ind_154animals_dPTA.Body_depth.txt | cut -f 2 | sed '1i#IID' > sample.list
ls ind_154animals_dPTA.*.txt | 
    while read id; do
    idd=$(echo $id | cut -d. -f2);
    cat $id | cut -f 3 | sed "1i$idd" > trait.tmp.$idd
    done
paste -d "\t" sample.list trait.tmp.* > cdcb.hol.traits
sed -i 's/\r//g' cdcb.hol.traits 

cp i154* cdcb.hol.traits /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/plink_gwas
cp cdcb.hol.traits /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/plink_gwas/ph.hol.traits
### # IID ph.hol.traits

cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/plink_gwas

plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-sv.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.01 --make-bed --out All2_BTAN_WGS.sv.bfile.maf1

plink2  --threads $nthreads --vcf ../3.pan_all/All2_BTAN_WGS-var.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --mind 0.1 --maf 0.01 --make-bed --out All2_BTAN_WGS.var.bfile.maf1

plink2  --threads $nthreads --vcf ../7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    --max-alleles 2 --const-fid 0 --chr-set 29 --chr 1-29 \
    --keep-allele-order --maf 0.01 --make-bed --out All2_BTAN_WGS.gsnp.bfile.maf1

plink2  --threads $nthreads --bfile All2_BTAN_WGS.gsnp.bfile.maf1 \
    --const-fid 0 --chr-set 29 \
    --pca --out All2_BTAN_WGS.gsnp.bfile.maf1

plink2  --threads $nthreads --bfile  All2_BTAN_WGS.var.bfile.maf1 \
    --const-fid 0 --chr-set 29 \
    --pca --out  All2_BTAN_WGS.var.bfile.maf1


###
cat > plink_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6

chrom=29
nthreads=$SLURM_CPUS_PER_TASK

if [[ $cov == "F" ]]; then
    plink2 --threads $nthreads --bfile $pre_gdata \
        --chr-set $chrom \
        --pheno iid-only $ph --pheno-name $trait \
        --no-categorical --glm allow-no-covars --out $out_pre.pc$pc
else
    plink2 --threads $nthreads --bfile $pre_gdata --chr-set 29 \
        --pheno iid-only $ph --pheno-name $trait \
        --covar $cov --covar-name $(printf "PC%d " $(seq 1 $pc)) \
        --no-categorical --glm --out $out_pre.pc$pc
fi
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

nthreads=16
for var in sv var gsnp; do

pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
out_pre=cdcb_gwas_nocv.$var
trait=Milk
cov=F
pc=F
ph=ph.hol.traits
#cov=All2_BTAN_WGS.var.bfile.maf1.eigenvec
#pc="2"

for trait in fat milk pro; do

sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash plink_gwas.sh $pre_gdata $out_pre $trait $cov $pc $ph
"
done
done
