#!/bin/bash

set -eu

panref_prefix=/90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes/ref_Holre_pangenie-var/ref_output
snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/svwgssnps/snp_input.$chr.maf01.vcf.gz

work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs
samples_chunk=$work_dir/tmp/chunk/sample_chunk_
out_prefix=$work_dir/umd50ks_HolrePan

chr_set=29
SBATCH_ACCOUNT=${SBATCH_ACCOUNT:-"bull_scr"}
nthreads=$SLURM_CPUS_PER_TASK

beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

mkdir -p $work_dir/logs $work_dir/tmp/chunk $panref_path/tmp
cd $work_dir

index=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
sample_index=$samples_chunk$index

# csplit -n 3 -s -f tmp/chunk/sample_chunk_ sample.list 100 {502}

main () {

    ref_chrSplit
    ref_beagle
    ref_rm
    input_chrSplit
    chrBeagle
    input_chunk
    chunkBeagle
    pmerge_chr
    bmerge_chunk
}

# Submit all jobs and wait for them to finish
job_ids=()

# Submit reference jobs
for chr in $(seq 1 $chr_set); do
    job_id=$(sbatch_sub -f "run_refChr" -c "$chr" -t 16)
    job_ids+=($job_id)
done

# Wait for reference jobs to finish
wait

# Submit imputation jobs
for chr in $(seq 1 $chr_set); do
    job_id=$(sbatch_sub -f "run_impChr" -c "$chr" -t 4)
    job_ids+=($job_id)
done

# Wait for imputation jobs to finish
wait

# Submit chunk jobs
for chr in $(seq 1 $chr_set); do
    job_id=$(sbatch_sub -f "run_impChunk" -c "$chr" -t 4 -a 502)
    job_ids+=($job_id)
done

# Wait for chunk jobs to finish
wait

# Submit merge jobs
for chr in $(seq 1 $chr_set); do
    job_id=$(sbatch_sub -f "run_merge" -c "$chr" -t 4)
    job_ids+=($job_id)
done

# Wait for merge jobs to finish
wait

sbatch_sub() {
    local func=""
    local nt=""
    local nc=""
    local dep=""
    local array=""
    local sbatch_options=()

    while getopts ":f:c:t:d:a:" opt; do
        case $opt in
            f)
                func="$OPTARG"
                ;;
            c)
                nc="$OPTARG"
                ;;
            t)
                nt="$OPTARG"
                sbatch_options+=("--cpus-per-task=$OPTARG")
                ;;
            d)
                dep="$OPTARG"
                sbatch_options+=("--dependency=afterok:$OPTARG")
                ;;
            a)
                array="$OPTARG"
                sbatch_options+=("--array=0-$OPTARG")
                ;;
            \?)
                echo "Invalid option: -$OPTARG" >&2
                return 1
                ;;
            :)
                echo "Option -$OPTARG requires an argument." >&2
                return 1
                ;;
        esac
    done

    shift $((OPTIND - 1))

    if [[ -z "$func" ]]; then
        echo "Error: Function name (-f) is required." >&2
        return 1
    fi

    mkdir -p logs

    local log_err="logs/$func.$nc.err"
    local log_out="logs/$func.$nc.out"

    # Construct the sbatch command using an array
    local sbatch_command=("sbatch")
    sbatch_command+=("--job-name=$func")
    sbatch_command+=("--chdir=$PWD")
    sbatch_command+=("${sbatch_options[@]}")
    sbatch_command+=("--error=$log_err")
    sbatch_command+=("--output=$log_out")

    # Construct the wrap command carefully
    local wrap_command="source config; $func"
    sbatch_command+=("--wrap=$wrap_command")

    if [[ "$DEBUG" == "1" ]]; then
        echo "${sbatch_command[@]}"
    fi

    # Execute sbatch and capture the job ID
    local job_id=$("${sbatch_command[@]}" | awk '{print $4}')

    if [[ -n "$job_id" ]]; then
        echo "Submitted job ID: $job_id"
        echo $job_id
    else
        echo "Error submitting job." >&2
        return 1
    fi

    return 0
}

function run_ref () {

    ref_beagle
}

function run_refChr () {
    ref_chrSplit
    ref_beagleChr
}

function run_impChr () {
    input_chrSplit
    chrBeagle
}

function run_impChunk () {
    input_chunk
    chunkBeagle
}

function run_merge () {
    bmerge_chunk
}

function run_gather () {
    pmerge_chr
}

function ref_beagle () {
    java -jar $beagle gt=$panref_vcf \
        out=$panref_prefix &&
    tabix -f -p vcf $panref_prefix.vcf.gz
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

function input_chunk () {
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

function bmerge_chunk () {
    bcftools merge --threads $nthreads \
        tmp/chunk/sample_*.afterImput.out.$chr.maf01.vcf.gz \
        -Oz -o tmp/sample.afterImput.out.$chr.vcf.gz &&
        tabix -f -p vcf tmp/sample.afterImput.out.$chr.vcf.gz
}

function pmerge_chr () {
    seq 1 $chr_set | while read chr; do
        echo tmp/sample.afterImput.out.$chr
    done > merge.list

    seq 1 $chr_set | while read chr; do
        plink2 --threads $nthreads tmp/sample.afterImput.out.$chr.vcf.gz \
            --chr-set $chr_set --const-fid --make-pgen --out tmp/sample.afterImput.out.$chr
    done

    plink2 --threads $nthreads --pmerge-list merge.list \
        --chr-set $chr_set --const-fid \
        --make-pgen --out $out_prefix

    plink2 --threads $nthreads --pfile $out_prefix --mind 0.1 --geno 0.05 --maf 0.01 \
        --chr-set $chr_set --const-fid --make-pgen --out $out_prefix.filter

}


