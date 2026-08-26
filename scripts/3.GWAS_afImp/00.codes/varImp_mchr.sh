#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

set -eu


chr_set=29
nthreads=$SLURM_CPUS_PER_TASK

work_dir=${PROJECT_ROOT}/hd_gwas/genotypes/

cd $work_dir

seq 1 $chr_set | while read chr; do
    echo sample.afterImput.out.$chr
done > merge.list

seq 1 $chr_set | while read chr; do
    plink2 --threads $nthreads sample.afterImput.out.$chr.vcf.gz \
        --chr-set $chr_set --const-fid --make-pgen --out sample.afterImput.out.$chr
done

plink2 --threads $nthreads --pmerge-list merge.list \
    --chr-set $chr_set --const-fid \
    --make-pgen --out sample.afterImput.out

plink2 --threads $nthreads --pfile sample.afterImput.out --mind 0.1 --geno 0.05 --maf 0.01 \
    --chr-set $chr_set --const-fid --make-pgen --out sample.afterImput.filter

rm sample.afterImput.out.*.vcf.gz*