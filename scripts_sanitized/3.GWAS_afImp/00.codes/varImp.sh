#!/bin/bash

set -eu

chr=$1
nthreads=$SLURM_CPUS_PER_TASK

snp_vcf=/90daydata/bull_age/liu.yang/data2024/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.$chr.vcf.gz
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/
panref_path=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/ref_Holre_pangenie-var/

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar


mkdir -p $work_dir
cd $work_dir

bcftools view -q 0.01:minor -r $chr --threads $nthreads \
    $snp_vcf \
    -oZ -o snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf snp_input.$chr.maf01.vcf.gz

java -jar $beagle ref=$panref_path/ref_output.$chr.vcf.gz \
    gt=snp_input.$chr.maf01.vcf.gz \
    out=snp_input.afterImput.out.$chr.maf01 &&
    tabix -f -p vcf snp_input.afterImput.out.$chr.maf01.vcf.gz

rm snp_input.$chr.maf01.vcf.gz*

sleep 0.2
