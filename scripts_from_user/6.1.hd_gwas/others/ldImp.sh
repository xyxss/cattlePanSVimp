
cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes

#cat > ldImp.sh << "EOF"
#!/bin/bash

set -eu

chr=$1
ref_type=$2
nthreads=$SLURM_CPUS_PER_TASK

work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/$ref_type/
if [ $ref_type == "svwgssnps" ]; then
    snp_vcf=/90daydata/bull_age/liu.yang/data2024/data_yahui_2024/umd/imputation/target_pop/hol_bull_imputation/allseq_1kbulls.hol.$chr.vcf.gz
else
    snp_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.maf01.reheader.vcf.gz
fi

panref_vcf=/90daydata/bull_age/liu.yang/imputation/hd_gwas/sv_gwas/Holsteinandrelated.pangenie-$ref_type.vcf.gz
haploid2diploid=/90daydata/bull_age/liu.yang/imputation/imputation/00.codes/haploid2diploid.awk
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

samples_chunk=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/samples/sample_chunk_

mkdir -p $work_dir
cd $work_dir

bcftools view -q 0.01:minor -r $chr --threads $nthreads \
    $snp_vcf \
    -oZ -o snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf snp_input.$chr.maf01.vcf.gz

bcftools view --threads $nthreads -r $chr \
    $panref_vcf | awk -f $haploid2diploid |
    bgzip -c > ref_input.$chr.vcf.gz &&
    tabix -f -p vcf ref_input.$chr.vcf.gz

java -jar $beagle gt=ref_input.$chr.vcf.gz \
    out=ref_output.$chr &&
    tabix -f -p vcf ref_output.$chr.vcf.gz

sbatch -A bull_scr \
    --array=0-502 \
    --job-name=$chr.${ref_type}Array \
    --cpus-per-task=16 \
    --error=logs/$ref_type.chr$chr.%a.err \
    --output=logs/$ref_type.chr$chr.%a.out \
    --wrap "
    ../imput.array.sbatch.sh $chr $work_dir
    "
    sleep 0.2
EOF

## imput.array.sbatch.sh

cat > imput.array.sbatch.sh << 'EOF2'
#!/bin/bash

set -eu

chr=$1
work_dir=$2
nthreads=$SLURM_CPUS_PER_TASK
samples_chunk=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/samples/sample_chunk_
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar


mkdir -p $work_dir/chunk/
cd $work_dir

index=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})

sample_index=$samples_chunk$index
bcftools view --threads $nthreads -r $chr \
    -S $sample_index --force-samples \
    snp_input.$chr.maf01.vcf.gz \
    -oZ -o chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz &&
    tabix -f -p vcf chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz

java -jar $beagle ref=ref_output.$chr.vcf.gz \
    gt=chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz \
    out=sample_${index}.afterImput.out.$chr.maf01 &&
    tabix -f -p vcf sample_${index}.afterImput.out.$chr.maf01.vcf.gz

rm  chunk/sample_${index}.snp_input.$chr.maf01.vcf.gz*

EOF2

## merge.array.sbatch.sh

cat > merge.array.sbatch.sh << 'EOF3'
#!/bin/bash

set -eu

work_dir=$2
nthreads=$SLURM_CPUS_PER_TASK
samples_chunk=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/samples/sample_chunk_
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar

cd $work_dir

bcftools concat --threads $nthreads \
    --allow-overlaps \
    sample_*.afterImput.out.$chr.maf01 
    -Oz -o sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf sample.afterImput.out.$chr.vcf.gz

EOF3



cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes

ref_type=svldsnps
ref_type=svwgssnps

seq 25 29 | while read chr; do
    sbatch -A bull_age \
        --job-name=chr$chr.$ref_type \
        --cpus-per-task=1 \
        --error=logs/$ref_type.chr$chr.err \
        --output=logs/$ref_type.chr$chr.out \
        --wrap="
        bash ldImp.sh $chr $ref_type
        "
done

ref_type=svldsnps
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/$ref_type/
echo 23 | while read chr; do
    sbatch \
         -A bull_scr \
        --array=0-502 \
        --job-name=$chr.${ref_type}Array \
        --cpus-per-task=8 \
        --error=logs/$ref_type.chr$chr.%a.err \
        --output=logs/$ref_type.chr$chr.%a.out \
        --wrap "
        bash imput.array.sbatch.sh $chr $work_dir
        "
done

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes
ref_type=svwgssnps
work_dir=/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/$ref_type/
seq 24 29 | while read chr; do
ind=$(sq | grep chr$chr.$ref_type | awk '{print $1}')
    sbatch -A bull_age --dependency=afterok:$ind \
         -A bull_scr \
        --array=0-502 \
        --job-name=$chr.${ref_type}Array \
        --cpus-per-task=8 \
        --error=logs/$ref_type.chr$chr.%a.err \
        --output=logs/$ref_type.chr$chr.%a.out \
        --wrap "
        bash imput.array.sbatch.sh $chr $work_dir
        "
done

ls svldsnps/sample_*.afterImput.out.*.maf01.vcf.gz | 
    cut -d. -f4 | sed -e 's|.*/||' | 
    sort | uniq -c  | tee svldsnps.check

ls svwgssnps/sample_*.afterImput.out.*.maf01.vcf.gz | 
    cut -d. -f4 | sed -e 's|.*/||' | 
    sort | uniq -c | tee svwgssnps.check

ls svldsnps/snp_input.*.maf01.vcf.gz | wc -l

ls svwgssnps/sample.*.vcf.gz

sbatch -A bull_scr \
    --array=0-502 \
    --job-name=$chr.${ref_type}Array \
    --cpus-per-task=16 \
    --error=logs/$ref_type.chr$chr.%a.err \
    --output=logs/$ref_type.chr$chr.%a.out \
    imput.array.sbatch.sh $chr $work_dir
    sleep 0.2
EOF

bcftools concat --threads $nthreads \
    --allow-overlaps \
    sample_*.afterImput.out.$chr.maf01.vcf.gz \
    -Oz -o sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf sample.afterImput.out.$chr.vcf.gz

nthreads=4
seq 1 29| while read chr; do
    sbatch -A bull_scr \
        --job-name=chr$chr.merge \
        --cpus-per-task=4 \
        --error=merge.chr$chr.err \
        --output=merge.chr$chr.out \
        --wrap="
bcftools merge --threads $nthreads \
    sample_*.afterImput.out.$chr.maf01.vcf.gz \
    -Oz -o sample.afterImput.out.$chr.vcf.gz &&
    tabix -f -p vcf sample.afterImput.out.$chr.vcf.gz
        "
done

ind=16576615
sbatch -A bull_scr --dependency=afterok:$ind  \
        --job-name=ch.merge \
        --cpus-per-task=16 \
        --error=merge.chr.err \
        --output=merge.chr.out \
        --wrap="
bcftools concat --threads $nthreads \
    --allow-overlaps svwgssnps/sample.*.vcf.gz \
    -Oz -o umdsample.afterImput.out.vcf.gz &&
    tabix -f -p vcf umdsample.afterImput.out.vcf.gz

plink2 --threads 8 --vcf umdsample.afterImput.out.vcf.gz \
    --chr-set 29 --const-fid --make-pgen --out umdsample.afterImput.out
plink2 --pfile umdsample.afterImput.out --mind 0.1 --geno 0.05 --maf 0.01 \
    --chr-set 29 --const-fid --make-pgen --out umdsample.afterImput.filter

        "

seq 1 29 | while read chr; do
    sbatch -A bull_scr \
        --job-name=ch$chr.merge \
        --cpus-per-task=16 \
        --error=merge.chr$chr.err \
        --output=merge.chr$chr.out \
        --wrap="
plink2 --threads 16 --vcf sample.afterImput.out.$chr.vcf.gz \
    --chr-set 29 --const-fid --make-pgen --out sample.afterImput.out.$chr.p
plink2 --pfile sample.afterImput.out.$chr.p --mind 0.1 --geno 0.05 --maf 0.01 \
    --chr-set 29 --const-fid --make-pgen --out sample.afterImput.out.$chr.pfliter
        "
done

seq 1 29 | while read chr; do 
echo sample.afterImput.out.$chr.pfliter >> file_list.txt
echo sample.afterImput.out.$chr.p >> file_list2.txt
done

sbatch -A bull_scr \
        --job-name=ch.merge \
        --cpus-per-task=16 \
        --error=merge.chr.err \
        --output=merge.chr.out \
        --wrap="
plink2  --threads 16 --pmerge-list file_list.txt --chr-set 29 --const-fid \
    --make-pgen --out ../umdsample.svwgssnp.afterImput.filter
"



sbatch -A bull_scr \
        --job-name=ch.merge \
        --cpus-per-task=16 \
        --error=merge.chr.err \
        --output=merge.chr.out \
        --wrap="
plink2 --threads 16 --pmerge-list file_list2.txt --chr-set 29 --const-fid \
--make-pgen --out ../umdsample.svwgssnp.afterImput
"
