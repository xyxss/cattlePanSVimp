
ref_path=/90daydata/bull_age/liu.yang/ref
ref_fa=$ref_path/ARS_UCD_v2.0.fa
ref_rm=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat


#### data preparation
work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp
mkdir -p $work_dir && cd $work_dir

run_file=$work_dir/run_files
code=$work_dir/00.codes
mkdir -p $run_file && mkdir -p $code

nthreads=$SLURM_CPUS_PER_TASK

for chr in $(seq 1 29); do
    bcftools view --threads $nthreads holPub.pangenie-sv.filter.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    awk -f $code/addID.awk PG="PG" |
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz
done

nt=8
proc=imp_all.sh
for chr in $(seq 1 29); do
    comp="chr.Holstein750.pangenie-sv"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz "rate:X,tool:F,rm:F,gf:F,typ:F,gt:F"
"
done

chr=10

### for comparison pangenie-var pangenie-sv deepv-snps pangenie-snps

bcftools view --threads $nthreads holPub.pangenie-sv.filter.vcf.gz \
    -r $chr -S Holstein750 --force-samples -Oz -o $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-sv.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub.pangenie-var.filter.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-var.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-var.chr$chr.vcf.gz

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-sv.chr$chr.vcf.gz \
    holPub.deepv-snps.chr$chr.vcf.gz \
    -Oz -o holPub_Thin0k.deepv-snps.chr$chr.vcf.gz &&
    tabix -f -p vcf holPub_Thin0k.deepv-snps.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub_Thin0k.deepv-snps.chr$chr.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.deepv-snps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.deepv-snps.chr$chr.vcf.gz

bcftools view --threads $nthreads -v snps holPub.pangenie-var.chr$chr.vcf.gz \
    -Oz -o holPub.pangenie-snps.chr$chr.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-snps.chr$chr.vcf.gz

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-sv.chr$chr.vcf.gz \
    holPub.pangenie-snps.chr$chr.vcf.gz \
    -Oz -o holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz &&
    tabix -f -p vcf holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz

nt=8
proc=imp_all.sh
for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
    comp="all.Holstein750.$sub"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.$sub.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done


## for population Holstein250 HolsteinRelated250 Jersey250 Holstein750 MultiBreed750

function grp_vcf () {
    local vcf_file=$1
    local out_file=$2
    for grp in Holstein250 HolsteinRelated250 Jersey250 MultiBreed750; do
        bcftools view --threads $nthreads -S $grp $vcf_file \
            -Oz -o ${out_file/holPub/$grp} && 
        tabix -f -p vcf ${out_file/holPub/$grp}
        echo ${out_file/holPub/$grp}
    done
}

grp_vcf holPub.pangenie-sv.chr$chr.vcf.gz run_files/holPub_Thin0k.pangenie-sv.chr$chr.vcf.gz
grp_vcf holPub.pangenie-var.chr$chr.vcf.gz run_files/holPub_Thin0k.pangenie-var.chr$chr.vcf.gz
grp_vcf holPub_Thin0k.deepv-snps.chr$chr.vcf.gz run_files/holPub_Thin0k.deepv-snps.chr$chr.vcf.gz
grp_vcf holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz run_files/holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz


nt=8
proc=imp_all.sh
for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
for grp in Holstein250 HolsteinRelated250 Jersey250 MultiBreed750; do
    comp="grp.$grp.$sub"
    sbatch -A bull_scr \
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

## for MAF 0.01 0.05 0.01 0.1 0.2

function MAF_vcf () {
    local vcf_file=$1
    for MAF in 0.05 0.01 0.1 0.2; do
        bcftools view --threads $nthreads -i "F_MISSING < 0.1 && MAF > $MAF" $vcf_file \
            -Oz -o ${vcf_file/Holstein750_Thin0k/Holstein750_Thin0k_${MAF/0./0}}
        tabix -f -p vcf ${vcf_file/Holstein750_Thin0k/Holstein750_Thin0k_${MAF/0./0}}
        echo ${vcf_file/Holstein750_Thin0k/Holstein750_Thin0k_${MAF/0./0}}
    done
}
for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
    MAF_vcf run_files/Holstein750_Thin0k.${sub}.chr$chr.vcf.gz
    for MAF in 0.05 0.01 0.1 0.2; do
        comp="MAF.$MAF.$sub"
        sbatch -A bull_scr \
            -D $PWD \
            --export=ALL \
            -J $comp.chr$chr.${proc//.*} \
            -c $nt \
            -o logs/$comp.chr$chr.${proc//.*}.out \
            --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k_${MAF/0./0}.$sub.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done
done


## for type

function typ_vcf () {
    local vcf_file=$1
    for typ in DEL INS COMPLEX; do
        zcat $vcf_file | awk '$1 ~ "#"{print;next} {split($3,a,/:|-|_/);if(a[4] == typ || a[4] ~ "s"){print}}' typ=$typ | bgzip -c > ${vcf_file/Holstein750_Thin0k/Holstein750_Thin0k_$typ} &&
        tabix -f -p vcf ${vcf_file/Holstein750_Thin0k/Holstein750_Thin0k_$typ}
    done
}

for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
    typ_vcf run_files/Holstein750_Thin0k.${sub}.chr$chr.vcf.gz
    for typ in DEL INS COMPLEX; do
        comp="typ.$typ.$sub"
        sbatch -A bull_scr \
            -D $PWD \
            --export=ALL \
            -J $comp.chr$chr.${proc//.*} \
            -c $nt \
            -o logs/$comp.chr$chr.${proc//.*}.out \
            --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k_$typ.$sub.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done
done


## for thin deepv-snps 
function thin_var () {
    local in_vcf=$1
    local thin=$2
    vcftools --gzvcf $in_vcf --thin ${thin} --recode --recode-INFO-all --stdout | 
    bgzip -c > ${in_vcf/Thin0k/Thin$((thin/1000))k} &&
    tabix -f -p vcf ${in_vcf/Thin0k/Thin$((thin/1000))k}
}

for thin in 1000 5000 10000 50000; do
    thin_var $run_file/Holstein750_Thin0k.deepv-snps.chr$chr.vcf.gz $thin &
done
wait

nt=8
proc=imp_all.sh
for thin in 1 5 10 50; do
    comp="thin.Holstein750.deepv-snps.Thin${thin}k"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin${thin}k.deepv-snps.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done

for thin in 1000 5000 10000 50000; do
    thin_var $run_file/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz $thin &
done
wait

nt=8
proc=imp_all.sh
for thin in 1 5 10 50; do
    comp="thin.Holstein750.pangenie-snps.Thin${thin}k"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin${thin}k.pangenie-snps.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done


# for var snp indel snv


bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-sv.chr$chr.vcf.gz \
    holPub.deepv-indels.chr$chr.vcf.gz \
    -Oz -o holPub_Thin0k.deepv-indels.chr$chr.vcf.gz &&
    tabix -f -p vcf holPub_Thin0k.deepv-indels.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub_Thin0k.deepv-indels.chr$chr.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    awk -f $code/addID.awk PG="DV" |
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.deepv-indels.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.deepv-indels.chr$chr.vcf.gz


# deepvariant
bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-sv.chr$chr.vcf.gz \
    holPub.deepv.chr$chr.vcf.gz \
    -Oz -o holPub_Thin0k.deepv.chr$chr.vcf.gz &&
    tabix -f -p vcf holPub_Thin0k.deepv.chr$chr.vcf.gz


bcftools view --threads $nthreads holPub_Thin0k.deepv.chr$chr.vcf.gz \
    -r $chr -S Holstein750 --force-samples |
    awk -f $code/addID.awk PG="DV" |
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.deepv-snv.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.deepv-snv.chr$chr.vcf.gz


nt=8
proc=imp_all.sh
for var in snv indels; do
    comp="var.deepv-$var"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.deepv-$var.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done



##
zcat run_files/Holstein750_Thin0k.deepv-snps.chr$chr.vcf.gz | grep -v "#" | cut -f 1-2 > deepv_positions_list.txt

bcftools view --threads $nthreads -T deepv_positions_list.txt run_files/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz -Oz -o run_files/Holstein750_Thin0k.pangenie-deepv.chr$chr.vcf.gz && tabix -f -p vcf run_files/Holstein750_Thin0k.pangenie-deepv.chr10.vcf.gz

# 125678 pangenie-deepv
# 508654 pangenie-snps
# 190400 deepv-snps

zcat run_files/Holstein750_Thin0k.deepv-snps.chr$chr.vcf.gz | grep -v "#" | awk '{split($3,a,/:|-|_/);if(a[5] <= 50) {print $1"\t"$2}}'| cut -f 1-2 > deepv_snp_list.txt

bcftools view --threads $nthreads -T ^deepv_snp_list.txt run_files/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz -Oz -o run_files/Holstein750_Thin0k.pangenie-vdeepv.chr$chr.vcf.gz
# 382976 pangenie-vdeepv

nt=8
proc=imp_all.sh
for var in snps deepv vdeepv; do
    comp="var.pangenie-$var"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin0k.pangenie-$var.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done


## for thin deepv-snps 
function thin_var () {
    local in_vcf=$1
    local thin=$2
    vcftools --gzvcf $in_vcf --thin ${thin} --recode --recode-INFO-all --stdout | 
    bgzip -c > ${in_vcf/Thin0k/Thin$((thin/1000))k} &&
    tabix -f -p vcf ${in_vcf/Thin0k/Thin$((thin/1000))k}
}

for thin in 1000 5000 10000 50000; do
    thin_var $run_file/Holstein750_Thin0k.pangenie-snps.chr$chr.vcf.gz $thin &
done
wait


nt=8
proc=imp_all.sh
for thin in 1 5 10 50; do
    comp="thin.Holstein750.pangenie-snps.Thin${thin}k"
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $run_file/Holstein750_Thin${thin}k.pangenie-snps.chr$chr.vcf.gz rate:T,tool:T,rm:T,gf:T,typ:T,gt:T
"
done
