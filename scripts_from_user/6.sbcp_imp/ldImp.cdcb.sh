

cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval



cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.sbatch.sh << 'EOF2'
#!/bin/bash
set -eu

chr=$1
work_dir=$2
snp_vcf=$3
nthreads=$SLURM_CPUS_PER_TASK

panref_prefix=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/ref_Holre_pangenie-var/ref_output

chr_set=29
SBATCH_ACCOUNT=${SBATCH_ACCOUNT:-"dor_rna"}
nthreads=$SLURM_CPUS_PER_TASK

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

mkdir -p $work_dir/logs $work_dir/tmp
cd $work_dir


main () {
    run_impChunk
}

function run_impChunk () {
    input_filter
    chrBeagle
}


function input_filter () {
    bcftools view --threads $nthreads -r $chr \
        $snp_vcf |
        awk '
        BEGIN {FS=OFS="\t"} 
        $1 ~ /#/{print; next;} 
        {
            for(i=10;i<=NF;i++){
                if($i=="."){$i="./."}
            }
            print
        }' |
        bgzip -c > tmp/sample.snp_input.$chr.maf01.vcf.gz &&
        tabix -f -p vcf tmp/sample.snp_input.$chr.maf01.vcf.gz
}

function chrBeagle () {
    java -jar $beagle ref=$panref_prefix.$chr.vcf.gz \
        gt=tmp/sample.snp_input.$chr.maf01.vcf.gz \
        out=sample.afterImput.out.$chr.maf01 &&
        tabix -f -p vcf sample.afterImput.out.$chr.maf01.vcf.gz
}

main
EOF2


cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval/1.pan-snp

nthreads=32
work_dir=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval
snp_vcf=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz
seq 1 29 | while read chr; do
    sbatch \
        -A uvm_mckay \
        --job-name=$chr.cdcb \
        --cpus-per-task=$nthreads \
        --error=logs/cdcb.chr$chr.%a.err \
        --output=logs/cdcb.chr$chr.%a.out \
        --time=1-00:00:00 \
        --wrap "
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.sbatch.sh \
        $chr $work_dir $snp_vcf
        "
done



cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval

mkdir -p 0.pan-snp 1.gatk-snp 2.1kbull-snp 3.hd-snp 4.ld-snp 5.rna-snp

bcftools view --threads 8 -q 0.01:minor \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 0.pan-snp/All2_BTAN_WGS.pan-snp.vcf.gz &&
    tabix -f -p vcf 0.pan-snp/All2_BTAN_WGS.pan-snp.vcf.gz

awk 'FNR==NR { a[$1 FS $2]; next } ($1 FS $2) in a' \
    <(zcat /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz | grep -v "#") \
    <(zcat /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz | grep -v "#" | cut -f 1-2) \
    > 1.gatk-snp/holPub.filter_ind.gatk.snp.pos
# 11926255
# All_BTAN_WGS_gatk.snp.vcf.gz 21102718
# All2_BTAN_WGS-snps.vcf.gz 12408901 

bcftools view --threads 8 -q 0.01:minor \
    -R 1.gatk-snp/holPub.filter_ind.gatk.snp.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 1.gatk-snp/All2_BTAN_WGS.gatk-snp.vcf.gz &&
    tabix -f -p vcf 1.gatk-snp/All2_BTAN_WGS.gatk-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/allseq_1kbulls.hol.chrall.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 2.1kbull-snp/All2_BTAN_WGS.1kbull-snp.vcf.gz &&
    tabix -f -p vcf 2.1kbull-snp/All2_BTAN_WGS.1kbull-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/GCF_002263795.3_ARS-UCD2.0.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 3.hd-snp/All2_BTAN_WGS.hd-snp.vcf.gz &&
    tabix -f -p vcf 3.hd-snp/All2_BTAN_WGS.hd-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/hol.chip.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 4.ld-snp/All2_BTAN_WGS.ld-snp.vcf.gz &&
    tabix -f -p vcf 4.ld-snp/All2_BTAN_WGS.ld-snp.vcf.gz

awk 'FNR==NR { a[$1 FS $2]; next } ($1 FS $2) in a' <(zcat /90daydata/bull_age/liu.yang/data_save/rnaSeq_HiFisamples/rna_vcf/rna.hifi19.sample_chr.vcf.gz | grep -v "#") /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/allseq_1kbulls.hol.chrall.pos > holPub.filter_ind.g-rnasnps.pos

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/holPub.filter_ind.g-rnasnps.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-snps.vcf.gz \
    -Oz -o 5.rna-snp/All2_BTAN_WGS.rna-snp.vcf.gz &&
    tabix -f -p vcf 5.rna-snp/All2_BTAN_WGS.rna-snp.vcf.gz

i=0
nthreads=32
for var in pan-snp gatk-snp 1kbull-snp hd-snp ld-snp rna-snp; do
    work_dir=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval/$i.$var/
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval/$i.$var/All2_BTAN_WGS.$var.vcf.gz
    seq 1 29 | while read chr; do
        sbatch \
            -A uvm_mckay \
            --job-name=$chr.cdcb \
            --cpus-per-task=$nthreads \
            --error=logs/cdcb.chr$chr.%a.err \
            --output=logs/cdcb.chr$chr.%a.out \
            --time=1-00:00:00 \
            --wrap "
            bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.sbatch.sh \
            $chr $work_dir $snp_vcf
            "
    done
    i=$((i+1))
done

l */*sample.afterImput.out.*gz | awk '$5 ~ "M"{print $9}' | cut -d/ -f1 | sort | uniq -c

i=0
nthreads=8
for var in pan-snp gatk-snp 1kbull-snp hd-snp ld-snp rna-snp; do
sbatch -A bull_age \
    --job-name=mergeChr \
    --cpus-per-task=8 \
    --wrap="
    bcftools concat --threads $nthreads \
        --allow-overlaps $i.$var/sample.afterImput.out*vcf.gz | 
    bcftools view --threads $nthreads \
        -q 0.01:minor -Oz -o imp.pang.$i.$var.vcf.gz &&
    tabix -f -p vcf imp.pang.$i.$var.vcf.gz
    "
    i=$((i+1))
done


cd /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/2.imputation_gatk
mkdir -p 0.pan-snp 1.gatk-snp 2.1kbull-snp 3.hd-snp 4.ld-snp 5.rna-snp

bcftools view --threads 8 -q 0.01:minor \
    -R ../1.imputation_eval/1.gatk-snp/holPub.filter_ind.gatk.snp.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 0.pan-snp/All2_BTAN_WGS.pan-snp.vcf.gz &&
    tabix -f -p vcf 0.pan-snp/All2_BTAN_WGS.pan-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 1.gatk-snp/All2_BTAN_WGS.gatk-snp.vcf.gz &&
    tabix -f -p vcf 1.gatk-snp/All2_BTAN_WGS.gatk-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/allseq_1kbulls.hol.chrall.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 2.1kbull-snp/All2_BTAN_WGS.1kbull-snp.vcf.gz &&
    tabix -f -p vcf 2.1kbull-snp/All2_BTAN_WGS.1kbull-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/GCF_002263795.3_ARS-UCD2.0.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 3.hd-snp/All2_BTAN_WGS.hd-snp.vcf.gz &&
    tabix -f -p vcf 3.hd-snp/All2_BTAN_WGS.hd-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/hol.chip.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 4.ld-snp/All2_BTAN_WGS.ld-snp.vcf.gz &&
    tabix -f -p vcf 4.ld-snp/All2_BTAN_WGS.ld-snp.vcf.gz

bcftools view --threads 8 -q 0.01:minor \
    -R /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp/holPub.filter_ind.g-rnasnps.pos \
    /90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/7.gatk/All_BTAN_WGS_gatk.snp.vcf.gz \
    -Oz -o 5.rna-snp/All2_BTAN_WGS.rna-snp.vcf.gz &&
    tabix -f -p vcf 5.rna-snp/All2_BTAN_WGS.rna-snp.vcf.gz


i=0
nthreads=32
for var in pan-snp gatk-snp 1kbull-snp hd-snp ld-snp rna-snp; do

    work_dir=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/2.imputation_gatk/$i.$var/
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/2.imputation_gatk/$i.$var/All2_BTAN_WGS.$var.vcf.gz
    seq 1 29 | while read chr; do
        sbatch \
            -A uvm_mckay \
            --job-name=$chr.cdcb \
            --cpus-per-task=$nthreads \
            --error=logs/cdcb.chr$chr.%a.err \
            --output=logs/cdcb.chr$chr.%a.out \
            --time=1-00:00:00 \
            --wrap "
            bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.sbatch.sh \
            $chr $work_dir $snp_vcf
            "
    done
    i=$((i+1))
done

i=0
nthreads=8
for var in pan-snp gatk-snp 1kbull-snp hd-snp ld-snp rna-snp; do
sbatch -A bull_age \
    --job-name=mergeChr \
    --cpus-per-task=8 \
    --wrap="
    bcftools concat --threads $nthreads \
        --allow-overlaps $i.$var/sample.afterImput.out*vcf.gz | 
    bcftools view --threads $nthreads \
        -q 0.01:minor -Oz -o imp.pang.$i.$var.vcf.gz &&
    tabix -f -p vcf imp.pang.$i.$var.vcf.gz
    "
    i=$((i+1))
done



bcftools view --threads 8 -q 0.01:minor \
/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/All2_BTAN_WGS-sv.vcf.gz \
-Oz -o All2_BTAN_WGS.pang-sv.vcf.gz &&
tabix -f -p vcf All2_BTAN_WGS.pang-sv.vcf.gz

1.imputation_eval
ls */All2*gz | while read id; do sbatch -A bull_age --wrap="bcftools stats $id > $id.sta"; done
ls */All2*sta | while read id; do cat $id | grep 'number of SNPs:' |awk '{split(IDD,a,".");print a[3],$NF}' IDD=$id; done
pan-snp 11934283
gatk-snp 11492334
1kbull-snp 10632332
hd-snp 571555
ld-snp 67356
rna-snp 1539272
val.sv 82355


2.imputation_gatk
ls */All2*gz | while read id; do sbatch -A bull_age --wrap="bcftools stats $id > $id.sta"; done
ls */All2*sta | while read id; do cat $id | grep 'number of SNPs:' |awk '{split(IDD,a,".");print a[3],$NF}' IDD=$id; done
pan-snp 11492334
gatk-snp 16135920
1kbull-snp 12842651
hd-snp 606440
ld-snp 68557
rna-snp 1628389
R2 for SVs. SNP?

F1 for each sample, each type, length?



bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz "rate:X,tool:F,rm:F,gf:F,typ:F,gt:F"

