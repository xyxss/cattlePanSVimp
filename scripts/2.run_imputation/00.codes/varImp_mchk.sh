#!/bin/bash
# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
# --------------------------

set -eu

chr=$1

nthreads=$SLURM_CPUS_PER_TASK

work_dir=${PROJECT_ROOT}/hd_gwas/genotypes/

cd $work_dir

bcftools merge --threads $nthreads \
    sample_*.afterImput.out.$chr.maf01.vcf.gz \
    -Oz -o sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf sample.afterImput.out.$chr.vcf.gz
