# --- site configuration ---
# Copy config.sh.example to config.sh at the repository root, edit the paths,
# then `source config.sh` before running this script.
: "${PROJECT_ROOT:?PROJECT_ROOT is unset - see config.sh.example at the repository root}"
: "${SLURM_ACCOUNT:?SLURM_ACCOUNT is unset - see config.sh.example at the repository root}"
# --------------------------

cd ${PROJECT_ROOT}/imputation/4.hd_imp

awk ' FNR==NR{a[$4]=$1":"$2;next} {if(a[$2])$2=a[$2];print}' ARS-UCD2.0.bed ${PROJECT_ROOT}/hd_gwas/hol.chip.bim.bak > hol.chip.bim

awk '$2 ~ ":"{split($2,aa,":");print aa[1]"\t"aa[2]}' hol.chip.bim > hol.chip.pos
## 1922 id err, change it to ARS-UCD2.0
# 68706

## ldsnps
bcftools view --threads $nthreads -T hol.chip.pos \
    ../1.holPub_imp/holPub.pangenie-snps.vcf.gz |
    awk -f $code/haploid2diploid.awk |
    bgzip -c > holPub.pangenie-ldsnps.vcf.gz && 
tabix -f -p vcf holPub.pangenie-ldsnps.vcf.gz

## svldsnps
bcftools concat --threads $nthreads \
    --allow-overlaps \
    ../1.holPub_imp/holPub.pangenie-sv.vcf.gz \
    holPub.pangenie-ldsnps.vcf.gz \
    -Oz -o holPub.pangenie-svldsnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-svldsnps.vcf.gz


#plink ld
#xxxx
work_dir=${PROJECT_ROOT}/imputation/1.holPub_imp
run_file=$work_dir/run_files
code=$work_dir/00.codes
mkdir -p $run_file

mkdir -p $work_dir/imp_runs
chr=10

bcftools view --threads $nthreads holPub.pangenie-ldsnps.vcf.gz \
    -r $chr -S ../1.holPub_imp/Holstein750 --force-samples | 
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-ldsnps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-ldsnps.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub.pangenie-svldsnps.vcf.gz \
    -r $chr -S ../1.holPub_imp/Holstein750 --force-samples |
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-svldsnps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-svldsnps.chr$chr.vcf.gz

chr=10
nt=2
proc=imp_all.sh
for sub in pangenie-svldsnps pangenie-ldsnps; do
    comp="all.Holstein750.$sub"
    sbatch -A ${SLURM_ACCOUNT} \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.$sub.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done

## xxxxx

function grp_vcf () {
    local vcf_file=$1
    local out_file=$2
    for grp in Holstein250 HolsteinRelated250 Jersey250 MultiBreed750; do
        bcftools view --threads $nthreads -r $chr \
             -S ../1.holPub_imp/$grp $vcf_file |
             awk -f $code/haploid2diploid.awk |
            bgzip -c ${out_file/holPub/$grp} && 
        tabix -f -p vcf ${out_file/holPub/$grp}
        echo ${out_file/holPub/$grp}
    done
}

grp_vcf holPub.pangenie-ldsnps.vcf.gz run_files/holPub_Thin0k.pangenie-ldsnps.chr$chr.vcf.gz
grp_vcf holPub.pangenie-svldsnps.vcf.gz run_files/holPub_Thin0k.pangenie-svldsnps.chr$chr.vcf.gz


chr=10
nt=2
proc=imp_all.sh
for sub in pangenie-ldsnps pangenie-svldsnps; do
for grp in Holstein250 HolsteinRelated250 Jersey250 MultiBreed750; do
    comp="grp.$grp.$sub"
    sbatch -A ${SLURM_ACCOUNT} \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/${grp}_Thin0k.$sub.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done
done

