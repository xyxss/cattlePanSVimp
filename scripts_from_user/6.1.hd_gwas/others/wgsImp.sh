

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes

#cat > wgsImp.sh << "EOF"
#!/bin/bash

set -eu

chr=$1
nthreads=$SLURM_CPUS_PER_TASK

work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/wgsImp/
snp_vcf=/90daydata/bull_age/liu.yang/data2024/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.$chr.vcf.gz
panref_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/sv_gwas/Holsteinandrelated.pangenie-svwgssnps.vcf.gz
haploid2diploid=/90daydata/bull_age/liu.yang/imputation/imputation/00.codes/haploid2diploid.awk
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

samples_chunk=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/samples/sample_chunk_

mkdir -p $work_dir
cd $work_dir

bcftools view --threads $nthreads -q 0.01:minor -r $chr \
    $snp_vcf \
    -Oz -o snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf snp_input.$chr.maf01.vcf.gz

bcftools view --threads $nthreads -r $chr \
    $panref_vcf | awk -f $haploid2diploid |
    bgzip -c > ref_input.$chr.vcf.gz &&
    tabix -f -p vcf ref_input.$chr.vcf.gz

java -jar $beagle gt=ref_input.$chr.vcf.gz \
    out=ref_output.$chr

sbatch -A bull_scr \
    --array=0-502 \
    --job-name=$chr.wgsImp \
    --cpus-per-task=2 \
    --error=logs/wgsImp.chr$chr.%a.err \
    --output=logs/wgsImp.chr$chr.%a.out \
    imput.array.sbatch.sh $chr $work_dir
    sleep 0.2

    "
    sleep 0.2
done
done
EOF

seq 1 29 | while read chr; do
    sbatch -A bull_scr \
        --job-name=chr$chr.wgsImp \
        --cpus-per-task=48 \
        --error=logs/wgsImp.chr$chr.err \
        --output=logs/wgsImp.chr$chr.out \
        --wrap="
    bash wgsImp.sh $chr
    "
    sleep 0.1
done



seq 1 23 | while read chr; do
#jobid=$(sq | grep chr$chr.wgsImp | awk '{print $1}')
    sbatch --dependency=afterany:$jobid \
        -A bull_scr \
        --job-name=chr$chr.wgsImp \
        --cpus-per-task=48 \
        --error=logs/wgsImp.chr$chr.err \
        --output=logs/wgsImp.chr$chr.out \
        --wrap="
    bash wgsImp.sh $chr
    "
    sleep 0.1
done

sbatch --dependency=afternotok:12345 second_job.slurm