

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/1.svImpCompWGS

#### cdcb WGS GWAS for SV&SNP ==================================================================================================================================

## SNP dataset

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

#### cdcb 173 GWAS for SV&SNP ==================================================================================================================================


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/1.svImpCompWGS

# phenotype
sed 's/\t/,/g' /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/plink_gwas/ph.hol.traits > ph.hol.traits.csv

# data

awk 'NR==FNR{a[$1,$2,$4,$5]=$3;next} {if(a[$1,$2,$4,$5] != "") {print $3,"===",a[$1,$2,$4,$5]}}' \
    ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8 \
    <(zcat /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-var.vcf.gz) |
    grep -v "##" > cdcb173wgs2umd50kVar.IDcovert

awk 'BEGIN{FS=OFS="\t"} 
    NR==FNR{a[$1,$2,$4,$5]=$3;next} 
    {if($1 ~ "#" || a[$1,$2,$4,$5] != "") {$3=a[$1,$2,$4,$5];print $0}}' \
    ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8 \
    <(zcat /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-var.vcf.gz) |
    bgzip -c > cdcb173wgs-var.vcf.gz

plink2  --threads 8 --vcf cdcb173wgs-var.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 --extract ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8 \
    --mind 0.1 --maf 0.01 --make-pgen --out cdcb173WGS_HolrePan.filter.r8

plink2  --threads 8 --vcf cdcb173wgs-var.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 --extract ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8 \
    --mind 0.1 --maf 0.01 --make-bed --out cdcb173WGS_HolrePan.filter.r8

cat cdcb173WGS_HolrePan.filter.r8.pvar | grep  -v  "#" | cut -f 3 |  awk '{split($1,a,/-|_/);print a[2]}' | sort | uniq -c
   6222 COMPLEX
  15054 DEL
  26073 INS
10530233 SNP

# data

nthreads=42
pre_gdata=cdcb173WGS_HolrePan.filter.r8
#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cat $pre_gdata.pvar | awk 'BEGIN{print "SNP"} $0 !~ "#"{print $3}' > $pre_gdata.snp_info.csv

grm=F
cov=F
pc=F
ph=ph.hol.traits.csv

for trait in milk fat pro Stature Body_depth; do
out_pre=cdcb-173WGS_slemm_$trait

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
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt
done


####


#### cdcb 173 GWAS for SV only ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/1.svImpCompWGS

plink2  --threads 8 --vcf cdcb173wgs-var.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 --extract ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8.sv \
    --mind 0.1 --maf 0.01 --make-pgen --out cdcb173WGSsv_HolrePan.filter.r8

plink2  --threads 8 --vcf cdcb173wgs-var.vcf.gz \
    --const-fid 0 --chr-set 29 --chr 1-29 --extract ../0.svImp_r8/umd50ksWGS_HolrePan.filter.r8.sv \
    --mind 0.1 --maf 0.01 --make-bed --out cdcb173WGSsv_HolrePan.filter.r8


nthreads=42
pre_gdata=cdcb173WGSsv_HolrePan.filter.r8
#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cat $pre_gdata.pvar | awk 'BEGIN{print "SNP"} $0 !~ "#"{print $3}' > $pre_gdata.snp_info.csv

grm=F
cov=F
pc=F
ph=ph.hol.traits.csv

for trait in milk fat pro Stature Body_depth; do
out_pre=cdcb-173WGSsv_slemm_$trait

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


pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6

nthreads=$SLURM_CPUS_PER_TASK
slemm_path=/project/bull_scr/liu.yang/software/slemm-v0.90.1-x86_64-linux/

$slemm_path/slemm --num_threads $nthreads --bfile $pre_gdata \
    --max_heritability 0.4 \
    --phenotype_file $ph --trait $trait \
    --lmm --snp_info_file $pre_gdata.snp_info.csv \
    --min_maf 0.01 \
    --seed 9 \
    --output_file $out_pre.pc$pc

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt
done

