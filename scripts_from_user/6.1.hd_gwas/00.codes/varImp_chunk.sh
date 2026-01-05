#!/bin/bash

set -eu

chr=$1
nthreads=$SLURM_CPUS_PER_TASK

snp_vcf=/90daydata/bull_age/liu.yang/data2024/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.$chr.vcf.gz
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/
panref_path=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/ref_Holre_pangenie-var/

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

samples_chunk=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/samples/sample_chunk_

mkdir -p $work_dir/chunk/
cd $work_dir

index=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
sample_index=$samples_chunk$index

bcftools view --threads $nthreads -r $chr \
    -S $sample_index --force-samples \
    snp_input.$chr.maf01.vcf.gz \
    -oZ -o chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz

java -jar $beagle ref=$panref_path/ref_output.$chr.vcf.gz \
    gt=chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz \
    out=sample_${index}.afterImput.out.$chr.maf01 &&
    tabix -f -p vcf sample_${index}.afterImput.out.$chr.maf01.vcf.gz

rm  chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz*
