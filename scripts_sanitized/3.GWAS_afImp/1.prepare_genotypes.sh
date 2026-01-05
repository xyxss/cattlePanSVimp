
mkdir -p /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes


# WGS SNP yahui 
seq 1 29 | while read chr; do


ls /90daydata/bull_age/liu.yang/data2024/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.*.vcf.gz | while read id; do 
    idd=$(basename $id);
sbatch -A bull_scr \
        --cpus-per-task=4 \
        --mem-per-cpu=8g \
        --wrap="
    bcftools view -q 0.01:minor --threads $nthreads $id  ${idd/.vcf.gz/.maf01.vcf.gz}
    "
done

ls *bed | sed 's/.bed//' > merge.list

plink2 --threads 8 --pmerge-list merge.list bfile \
    --chr-set 29 --const-fid --maf 0.01 \
    --make-pgen --out allseq_1kbulls.hol.chrall.maf01

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/sv_gwas

### give reference for Hol and Hol-related samples
vcf=/90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/holPub.pangenie-svwgssnps.vcf.gz
cat /90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/holPub.samples.group | awk -F "\t" '$2 ~ "Holstein" {print $1}' > Holsteinandrelated.samples

bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svwgssnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svwgssnps.vcf.gz

vcf=/90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/holPub.pangenie-svldsnps.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svldsnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svldsnps.vcf.gz

vcf=/90daydata/bull_age/liu.yang/imputation/imputation/5.rna_shp/holPub.pangenie-svrnassnps.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svrnassnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svrnassnps.vcf.gz

!!!!!!!
vcf=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/holPub.pangenie-var.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-var.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-var.vcf.gz


plink2 --bfile hol.chip --chr-set 29 --const-fid --maf 0.01 -recode vcf --out hol.chip.maf01


zcat hol.chip.maf01.vcf.gz | awk '$1 ~ "#"{print;next}{exit}'| grep -v "chrSet" > cleaned_header.txt
bcftools view -H hol.chip.maf01.vcf.gz > data_no_header.vcf
bcftools reheader  --threads 8  -h cleaned_header.txt hol.chip.maf01.vcf.gz -oZ -o hol.chip.maf01.reheader.vcf.gz

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/
mkdir -p samples
csplit -n 3 -s -f samples/sample_chunk_ sample.list 100 {502}


 awk '{print $1,$1":"$4,$3,$4,$5,$6}' allseq_1kbulls.hol.chrall.maf01.bim.bak > allseq_1kbulls.hol.chrall.maf01.bim
 


#### 


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs
bfile=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs/umd50ksWGS_HolrePan.filter
cat $bfile.bim | 
       awk '$2 !~ "SNP" && $2 !~ "s"{print $2}' > $(basename $bfile).sv
plink2 --threads 8 --pfile $bfile --extract $(basename $bfile).sv --chr-set 29 --make-bed --out umd50ksWGSsv_HolrePan.filter


 cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/6.svImp_svld
bfile=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld/umd50ksLD_HolrePan.filter2
cat $bfile.bim | 
       awk '$2 !~ "SNP" && $2 !~ "s"{print $2}' > $(basename $bfile).sv
plink2 --threads 8 --bfile $bfile --extract $(basename $bfile).sv --chr-set 29 --make-bed --out umd50ksLDsv_HolrePan.filter


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes


awk '$7 ~ "IMP" && $3 ~ "SV"{split($7,a,/;|=/);if(a[3] > 0.85)print }' \
    umd50ksLD_EupSV.out.pvar > umd50ksLD_EupSV.out.filter

awk '$3 !~ "SNP" && $3 !~ "s"{print $2}' umdsample.svwgssnp.afterImput.filter.pvar

awk '$7 ~ "IMP" && $3 ~ "SV"{split($7,a,/;|=/);if(a[3] > 0.85)print }' \
    umdsample.svwgssnp.afterImput.filter.pvar > umdsample.svwgssnp.out.filter


cd /90daydata/bull_age/liu.yang/EuropeanCattleItems
awk '$3 ~ "SV"' umd50ksLD_EupSV.out.pvar > umd50ksLD_EupSV.out.sv
plink2 --threads 8 --bfile umd50ksLD_EupSV.filter --chr-set 29 --make-bed --out umd50ksLD_EupSV.sv
