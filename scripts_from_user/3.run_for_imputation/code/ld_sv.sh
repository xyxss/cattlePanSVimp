#!/bin/bash

#set -o nounset
#set -o errexit

comp=$1
chr=$2
input_vcf_gz=$3
nthreads=$SLURM_CPUS_PER_TASK
#pangenieSV
rtgNcpu=1


####
work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/ld_runs
code=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/00.codes

mkdir -p $work_dir/$comp/$comp.chr$chr
cd $work_dir/$comp/$comp.chr$chr


#plink2 --threads $nthreads --vcf $input_vcf_gz --freq --out ld.$chr.maf.out


plink --threads $nthreads --vcf $input_vcf_gz \
--const-fid 0 --r2 \
--ld-window-kb 1000 \
--ld-window-r2 0.05 \
--out ld.$chr \
--vcf-half-call m

awk -f $code/plink_ld_extaTypes.awk ld.$chr.ld > $work_dir/$comp/$comp.chr$chr.ld.$chr.ld.ext
awk '$9 !~ "snv-snv"'  $work_dir/$comp/$comp.chr$chr.ld.$chr.ld.ext > $work_dir/$comp/$comp.chr$chr.ld.$chr.ld.sv

