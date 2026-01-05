cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes


#cat > ldImp.sh << "EOF"
#!/bin/bash

set -eu

chr=$1

panref_prefix=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/ref_Holre_pangenie-var/ref_output
snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz

work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs
samples_chunk=$work_dir/tmp/chunk/sample_chunk_
out_prefix=$work_dir/umd50ks_HolrePan

chr_set=29
SBATCH_ACCOUNT=${SBATCH_ACCOUNT:-"dor_rna"}
nthreads=$SLURM_CPUS_PER_TASK

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

mkdir -p $work_dir/logs $work_dir/tmp/chunk
cd $work_dir

main () {
    run_refChr
    run_impChr
}

function run_refChr () {
    ref_chrSplit
    ref_beagleChr
}

function run_impChr () {
    input_chunk
    chunkBeagle
}

function ref_chrSplit () {
    bcftools view --threads $nthreads -r $chr -q 0.01:minor \
        $panref_vcf | 
        awk '
        BEGIN {FS=OFS="\t"} 
        $1 ~ /#/{print; next;} 
        {
            for(i=10;i<=NF;i++){
                if($i=="."){$i="./."}
            }
            print
        }' |
        bgzip -c > tmp/ref_input.$chr.vcf.gz &&
        tabix -f -p vcf tmp/ref_input.$chr.vcf.gz
}

function ref_beagleChr () {
    java -jar $beagle gt=tmp/ref_input.$chr.vcf.gz \
        out=$panref_prefix.$chr &&
        tabix -f -p vcf $panref_prefix.$chr.vcf.gz
    rm tmp/ref_input.$chr.vcf.gz*
}

function input_chrSplit () {
    bcftools view --threads $nthreads -r $chr -q 0.01:minor \
        $snp_vcf \
        -oZ -o tmp/snp_input.$chr.maf01.vcf.gz &&
        tabix -f -p vcf tmp/snp_input.$chr.maf01.vcf.gz
}

function chrBeagle () {
    java -jar $beagle ref=$panref_prefix.$chr.vcf.gz \
        gt=tmp/snp_input.$chr.maf01.vcf.gz \
        out=tmp/snp_input.afterImput.out.$chr.maf01 &&
        tabix -f -p vcf tmp/snp_input.afterImput.out.$chr.maf01.vcf.gz
}


main


EOF

seq 1 29 | while read chr; do
    sbatch \
         -A bull_age \
        --array=0-503 \
        --job-name=$chr.wgsArray \
        --cpus-per-task=8 \
        --error=logs/wgsArray.chr$chr.%a.err \
        --output=logs/wgsArray.chr$chr.%a.out \
        --wrap "
        bash ldImp.sh $chr
        "
done



## /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.array.sbatch.sh

cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.array.sbatch.sh << 'EOF2'
#!/bin/bash
set -eu

chr=$1
work_dir=$2
snp_vcf=$3
nthreads=$SLURM_CPUS_PER_TASK

panref_prefix=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/ref_Holre_pangenie-var/ref_output

samples_chunk=$work_dir/tmp/chunk/sample_chunk_
chr_set=29
SBATCH_ACCOUNT=${SBATCH_ACCOUNT:-"dor_rna"}
nthreads=$SLURM_CPUS_PER_TASK

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

mkdir -p $work_dir/logs $work_dir/tmp/chunk
cd $work_dir

index=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
sample_index=$samples_chunk$index

# csplit -n 3 -s -f tmp/chunk/sample_chunk_ sample.list 100 {502}

main () {
    run_impChunk
}

function run_impChunk () {
    input_chunk
    chunkBeagle
}


function input_chunk () {
    if [ ! -f $snp_vcf.tbi ]; then
        tabix -f -p vcf $snp_vcf
    fi
    bcftools view --threads $nthreads -r $chr \
        -S $sample_index --force-samples \
        $snp_vcf \
        -oZ -o tmp/chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz &&
        tabix -f -p vcf tmp/chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz
}

function chunkBeagle () {
    java -jar $beagle ref=$panref_prefix.$chr.vcf.gz \
        gt=tmp/chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz \
        out=tmp/chunk/sample_${index}.afterImput.out.$chr.maf01 &&
        tabix -f -p vcf tmp/chunk/sample_${index}.afterImput.out.$chr.maf01.vcf.gz
}

main
EOF2

cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.array.sbatch.sh << 'EOF3'
#!/bin/bash
set -eu
chr=$1
work_dir=$2
nthreads=$SLURM_CPUS_PER_TASK
chr_set=29
cd $work_dir
bcftools merge --threads $nthreads \
    tmp/chunk/sample_*.afterImput.out.$chr.maf01.vcf.gz \
    -Oz -o tmp/sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf tmp/sample.afterImput.out.$chr.vcf.gz

plink2 --threads $nthreads --vcf tmp/sample.afterImput.out.$chr.vcf.gz \
        --chr-set $chr_set --const-fid --make-pgen --out tmp/sample.afterImput.out.$chr

EOF3

cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.chr.sbatch.sh << 'EOF4'
#!/bin/bash
set -eu
out_prefix=$1
work_dir=$2
nthreads=$SLURM_CPUS_PER_TASK

chr_set=29
cd $work_dir
seq 1 $chr_set | while read chr; do
    echo tmp/sample.afterImput.out.$chr
done > merge.list

plink2 --threads $nthreads --pmerge-list merge.list \
    --chr-set $chr_set --const-fid \
    --make-bed --out $out_prefix.out

plink2 --threads $nthreads --bfile $out_prefix.out \
    --mind 0.1 --geno 0.05 --maf 0.01 \
    --chr-set $chr_set --const-fid \
    --make-bed --out $out_prefix.filter
EOF4


nthreads=8
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs
#snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz
seq 1 29 | while read chr; do
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz
    sbatch \
        -A uvm_mckay \
        --array=0-503 \
        --job-name=$chr.wgsArray \
        --cpus-per-task=$nthreads \
        --error=logs/wgsArray.chr$chr.%a.err \
        --output=logs/wgsArray.chr$chr.%a.out \
        --wrap "
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.array.sbatch.sh $chr $work_dir $snp_vcf
        " | awk '{print $4}' | xargs -I {} \
        sbatch -A uvm_mckay --dependency=afterok:{} \
        --job-name=chk$chr.merge \
        --cpus-per-task=$nthreads \
        --error=logs/merge.chk$chr.err \
        --output=logs/merge.chk$chr.out \
        --wrap="
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.array.sbatch.sh $chr $work_dir
        " 
        done 
        

    sbatch -A bull_age \
        --job-name=mergeChr \
        --cpus-per-task=$nthreads \
        --error=logs/mergeChr.err \
        --output=logs/mergeChr.out \
        --wrap="
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.chr.sbatch.sh \
        $work_dir/umd50ksWGS_HolrePan $work_dir
        "


## /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.array.sbatch.sh
####
nthreads=8
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/3.svImp_cCattle
mkdir -p $work_dir/tmp/chunk $work_dir/logs
#bcftools query -l /90daydata/bull_age/liu.yang/imputation/cGTEx_data/data/lfang2/Imputeded_SNPs/Chr1.all-beagle.vcf.gz  | wc -l
#csplit -n 3 -s -f tmp/chunk/sample_chunk_ sample.list 100 {72}

####
nthreads=8
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/3.svImp_cCattle
#snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz
seq 1 29 | while read chr; do
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/cGTEx_data/data/lfang2/Imputeded_SNPs/Chr$chr.all-beagle.vcf.gz
    sbatch \
        -A bull_scr \
        --array=0-72 \
        --job-name=$chr.cCArray \
        --cpus-per-task=$nthreads \
        --error=logs/cCArray.chr$chr.%a.err \
        --output=logs/cCArray.chr$chr.%a.out \
        --wrap "
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.array.sbatch.sh $chr $work_dir $snp_vcf
        " | awk '{print $4}' | xargs -I {} \
        sbatch -A bull_scr --dependency=afterok:{} \
        --job-name=chk$chr.cCArraymerge \
        --cpus-per-task=8 \
        --error=logs/cCArraymerge.chk$chr.err \
        --output=logs/cCArraymerge.chk$chr.out \
        --wrap="
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.array.sbatch.sh $chr $work_dir
        " 
    done 
        

    sbatch -A bull_age \
    --job-name=mergeChr \
    --cpus-per-task=8 \
    --error=logs/mergeChr.err \
    --output=logs/mergeChr.out \
    --wrap="
    bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.chr.sbatch.sh $work_dir/cattlegtex_HolrePan $work_dir
    "



work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld
snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.maf01.reheader.vcf.gz
mkdir -p $work_dir/tmp/chunk $work_dir/logs
cd $work_dir
bcftools query -l $snp_vcf > sample.list
csplit -n 3 -s -f tmp/chunk/sample_chunk_ sample.list 10000 {4}

nthreads=8
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld
#snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz
seq 13 29 | while read chr; do
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.maf01.reheader.vcf.gz
    sbatch \
        -A uvm_mckay \
        --array=0-4 \
        --job-name=$chr.ldArray \
        --cpus-per-task=$nthreads \
        --error=logs/ldArray.chr$chr.%a.err \
        --output=logs/ldArray.chr$chr.%a.out \
        --wrap "
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/imput.array.sbatch.sh $chr $work_dir $snp_vcf
        " | awk '{print $4}' | xargs -I {} \
        sbatch -A uvm_mckay --dependency=afterok:{} \
        --job-name=chk$chr.merge \
        --cpus-per-task=$nthreads \
        --error=logs/merge.chk$chr.err \
        --output=logs/merge.chk$chr.out \
        --wrap="
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.array.sbatch.sh $chr $work_dir
        " 
        done 
        

    sbatch -A bull_age \
        --job-name=mergeChr \
        --cpus-per-task=$nthreads \
        --error=logs/mergeChr.err \
        --output=logs/mergeChr.out \
        --wrap="
        bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/merge.chr.sbatch.sh $work_dir/umd50ksLD_HolrePan $work_dir
        "

out_prefix=cattlegtex_HolrePan
    sbatch -A bull_age \
        --job-name=mergeChr \
        --cpus-per-task=32 \
        --error=logs/mergeChr.err \
        --output=logs/mergeChr.out \
        --wrap="
plink2 --threads 32 --pfile $out_prefix.filter \
    --chr-set 29 --const-fid \
    --make-bed --out $out_prefix.filter2
"



