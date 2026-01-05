#!/bin/bash

set -eu

chr=$1
nthreads=$SLURM_CPUS_PER_TASK

work_dir=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/ref_Holre_pangenie-var
panref_vcf=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/Holsteinandrelated.pangenie-var.vcf.gz
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

mkdir -p $work_dir
cd $work_dir

bcftools view --threads $nthreads -r $chr \
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
    bgzip -c > ref_input.$chr.vcf.gz &&
    tabix -f -p vcf ref_input.$chr.vcf.gz

java -jar $beagle gt=ref_input.$chr.vcf.gz \
    out=ref_output.$chr &&
    tabix -f -p vcf ref_output.$chr.vcf.gz

sleep 0.2
rm ref_input.$chr.vcf.gz*
