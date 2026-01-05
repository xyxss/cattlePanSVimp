#!/bin/bash

set -eu

chr=$1
nthreads=$SLURM_CPUS_PER_TASK

snp_vcf=$2
work_dir=$3
panref_path=$4
samples_chunk=$5

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar


mkdir -p $work_dir/chunk/
cd $work_dir

index=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
sample_index=$samples_chunk$index

bcftools view --threads $nthreads -r $chr \
    -S $sample_index --force-samples \
    $snp_vcf \
    -oZ -o chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz

java -jar $beagle ref=$panref_path/ref_output.$chr.vcf.gz \
    gt=chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz \
    out=sample_${index}.afterImput.out.$chr.maf01 &&
    tabix -f -p vcf sample_${index}.afterImput.out.$chr.maf01.vcf.gz

rm  chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz*

