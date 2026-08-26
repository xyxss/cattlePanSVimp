# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${DATA_ROOT:?DATA_ROOT is unset - see config.sh.example at the repository root}"
: "${DATA_DIR:?DATA_DIR is unset - see config.sh.example at the repository root}"
: "${SAMPLE_GROUP_TABLE:?SAMPLE_GROUP_TABLE is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
# --------------------------

mkdir -p ${PROJECT_ROOT}/hd_gwas/genotypes
cd ${PROJECT_ROOT}/hd_gwas/genotypes


# WGS SNP yahui 


ls ${DATA_DIR}/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.*.vcf.gz | while read id; do 
    idd=$(basename $id);
sbatch -A ${SLURM_ACCOUNT} \
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

cd ${PROJECT_ROOT}/hd_gwas/sv_gwas

### give reference for Hol and Hol-related samples
vcf=${PROJECT_ROOT}/imputation/4.hd_imp/holPub.pangenie-svwgssnps.vcf.gz
cat ${SAMPLE_GROUP_TABLE} | awk -F "\t" '$2 ~ "Holstein" {print $1}' > Holsteinandrelated.samples

bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svwgssnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svwgssnps.vcf.gz

vcf=${PROJECT_ROOT}/imputation/4.hd_imp/holPub.pangenie-svldsnps.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svldsnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svldsnps.vcf.gz

vcf=${PROJECT_ROOT}/imputation/5.rna_shp/holPub.pangenie-svrnassnps.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-svrnassnps.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-svrnassnps.vcf.gz

!!!!!!!
vcf=${PROJECT_ROOT}/imputation/1.holPub_imp/holPub.pangenie-var.vcf.gz
bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-var.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-var.vcf.gz


plink2 --bfile hol.chip --chr-set 29 --const-fid --maf 0.01 -recode vcf --out hol.chip.maf01


zcat hol.chip.maf01.vcf.gz | awk '$1 ~ "#"{print;next}{exit}'| grep -v "chrSet" > cleaned_header.txt
bcftools view -H hol.chip.maf01.vcf.gz > data_no_header.vcf
bcftools reheader  --threads 8  -h cleaned_header.txt hol.chip.maf01.vcf.gz -oZ -o hol.chip.maf01.reheader.vcf.gz

cd ${PROJECT_ROOT}/hd_gwas/genotypes/
mkdir -p samples
csplit -n 3 -s -f samples/sample_chunk_ sample.list 100 {502}


 awk '{print $1,$1":"$4,$3,$4,$5,$6}' allseq_1kbulls.hol.chrall.maf01.bim.bak > allseq_1kbulls.hol.chrall.maf01.bim
 


#### 


cd ${PROJECT_ROOT}/hd_gwas/genotypes/5.svImp_svwgs
bfile=${PROJECT_ROOT}/hd_gwas/genotypes/1.svImp_wgs/umd50ksWGS_HolrePan.filter
cat $bfile.bim | 
       awk '$2 !~ "SNP" && $2 !~ "s"{print $2}' > $(basename $bfile).sv
plink2 --threads 8 --pfile $bfile --extract $(basename $bfile).sv --chr-set 29 --make-bed --out umd50ksWGSsv_HolrePan.filter


 cd ${PROJECT_ROOT}/hd_gwas/genotypes/6.svImp_svld
bfile=${PROJECT_ROOT}/hd_gwas/genotypes/2.svImp_ld/umd50ksLD_HolrePan.filter2
cat $bfile.bim | 
       awk '$2 !~ "SNP" && $2 !~ "s"{print $2}' > $(basename $bfile).sv
plink2 --threads 8 --bfile $bfile --extract $(basename $bfile).sv --chr-set 29 --make-bed --out umd50ksLDsv_HolrePan.filter


cd ${PROJECT_ROOT}/hd_gwas/genotypes


awk '$7 ~ "IMP" && $3 ~ "SV"{split($7,a,/;|=/);if(a[3] > 0.85)print }' \
    umd50ksLD_EupSV.out.pvar > umd50ksLD_EupSV.out.filter

awk '$3 !~ "SNP" && $3 !~ "s"{print $2}' umdsample.svwgssnp.afterImput.filter.pvar

awk '$7 ~ "IMP" && $3 ~ "SV"{split($7,a,/;|=/);if(a[3] > 0.85)print }' \
    umdsample.svwgssnp.afterImput.filter.pvar > umdsample.svwgssnp.out.filter


cd ${DATA_ROOT}/EuropeanCattleItems
awk '$3 ~ "SV"' umd50ksLD_EupSV.out.pvar > umd50ksLD_EupSV.out.sv
plink2 --threads 8 --bfile umd50ksLD_EupSV.filter --chr-set 29 --make-bed --out umd50ksLD_EupSV.sv
