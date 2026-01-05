#!/bin/bash

set -eu

chr=$1

nthreads=$SLURM_CPUS_PER_TASK

work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/

cd $work_dir

bcftools merge --threads $nthreads \
    sample_*.afterImput.out.$chr.maf01.vcf.gz \
    -Oz -o sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf sample.afterImput.out.$chr.vcf.gz
