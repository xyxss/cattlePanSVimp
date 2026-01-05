
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


function thin_var () {
    local vcf=$1
    local thin=$2
    local thin1k=$((thin/1000))
    vcftools --gzvcf $vcf --thin ${thin} --recode --recode-INFO-all --stdout | 
    bgzip -c > ${vcf/.vcf.gz/.thin${thin1k}k.vcf.gz} &&
    tabix -f -p vcf ${vcf/.vcf.gz/.thin${thin1k}k.vcf.gz}
}

function concat_vcf () {
    local in_vcf=$1
    local out_vcf=$2
    bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-sv.chr$chr.vcf.gz \
    $in_vcf \
    -Oz -o $out_vcf
}

function thin_concat () {
    local in_vcf=$1
    local thins=$2
    local fthin=$3
    local ncpu=$nthreads
    local count=0
    if [ "$fthin" == "T" ]; then
        for thin in $thins; do
            thin_var $in_vcf $thin
        done
    fi
    for thin in $thins; do
    concat_vcf ${in_vcf/.vcf.gz/.thin$((thin/1000))k.vcf.gz} \
        ${in_vcf/./_Thin$((thin/1000))k.}
    echo ${in_vcf/./_Thin$((thin/1000))k.}
    ((count++))
    if (( count % ncpu == 0 )); then
        wait
    fi
    done 
}

function proc_vcf () {
    local in_vcf=$1
    bcftools view --threads $nthreads -S Holstein $in_vcf \
            -Oz -o ${in_vcf/holPub/Holstein}
    grp_vcf $in_vcf
    MAF_vcf ${in_vcf/holPub/Holstein}
    typ_vcf ${in_vcf/holPub/Holstein}
    cp $in_vcf $in_vcf
    echo $in_vcf
}

function grp_vcf () {
    local vcf_file=$1
    for grp in Holstein Holstein250 HolsteinRelated250 Jersey250 Holstein750 MultiBreed750; do
        bcftools view --threads $nthreads -S $grp $vcf_file \
            -Oz -o ${vcf_file/holPub/$grp}
        echo ${vcf_file/holPub/$grp}
    done
}

function MAF_vcf () {
    local vcf_file=$1
    for MAF in 0.05 0.01 0.1 0.2; do
        bcftools view --threads $nthreads -i "F_MISSING < 0.1 && MAF > $MAF" $vcf_file \
            -Oz -o ${vcf_file/Holstein/Holstein.${MAF/0./0}}
        echo ${vcf_file/Holstein/Holstein.${MAF/0./0}}
    done
}

function typ_vcf () {
    local vcf_file=$1
    for typ in DEL INS MNP COMPLEX; do
        bcftools view --threads $nthreads -i "INFO/SVTYPE == \"$typ\"" $vcf_file \
            -Oz -o ${vcf_file/pangenie-sv/pangenie-$typ}
        echo ${vcf_file/pangenie-sv/pangenie-$typ}
    done
}

function Holstein_vcf () {
    local vcf_file=$1
    bcftools view --threads $nthreads -S Holstein $vcf_file \
        -Oz -o ${vcf_file/holPub/Holstein}
    echo ${vcf_file/holPub/Holstein}
}

# ================================================================================================ #
# ================================================================================================ #
# ================================================================================================ #
## match samples and chromosomes

cat holPub.samples.group | awk -F "\t" '{print $1}' > holPub.samples
cat holPub.samples.group | awk -F "\t" '$2 == "Holstein" {print $1}' > Holstein
cat holPub.samples.group | awk -F "\t" '$2 == "Holstein" {print $1}' | shuf -n 750 > Holstein750
#963
shuf -n 250 Holstein > Holstein250
cat holPub.samples.group | awk -F "\t" '$2 ~ "Holstein-X-Jersey" {print $1}' | shuf -n 250 > HolsteinRelated250
cat holPub.samples.group | awk -F "\t" '$2 == "Jersey" {print $1}' | shuf -n 250 > Jersey250
cat Holstein250 HolsteinRelated250 Jersey250 > MultiBreed750

# pangenie with all chromosomes
vcf=/90daydata/bull_age/liu.yang/imputation/pangenie_holPub/3.pan_all/hol-pg2hic-2024-05-22_graph_genotyping.merge-biallelic.vcf.gz

bcftools view --threads $nthreads -S holPub.samples \
    -i 'F_MISSING < 0.1 && AC > 0' $vcf |
    bcftools annotate --threads $nthreads --remove QUAL,INFO,FORMAT \
    -Oz -o holPub.pangenie.vcf.gz &&
    tabix -f -p vcf holPub.pangenie.vcf.gz

bcftools view --threads $nthreads -S holPub.samples \
    -i 'F_MISSING < 0.1 && MAF > 0.01' $vcf |
    bcftools annotate --threads $nthreads --remove QUAL,INFO,FORMAT \
    -Oz -o holPub.pangenie.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie.filter.vcf.gz

# deepvariant with all chromosomes

vcf2=/90daydata/bull_age/liu.yang/imputation/pangenie_holPub/4.deepv_all/ARS_UCD_v2.0.chr.fa.deepv-bi.vcf.gz
bcftools view --threads $nthreads -S holPub.samples \
    -i 'F_MISSING < 0.1 && MAF > 0.01' $vcf2 |
    bcftools annotate --threads $nthreads --remove QUAL,INFO,FORMAT |
    awk -f $code/addID.awk PG="DV" |
    bgzip -c -@ $nthreads > holPub.deepv.filter.vcf.gz && 
    tabix -f -p holPub.deepv.filter.vcf.gz

# pangenie-sv with all chromosomes
bcftools view --threads 4 holPub.pangenie.vcf.gz \
    -r $(seq 1 29 | tr '\n' ',') \
    -Oz -o holPub.pangenie.vcf2.gz && 
    tabix -f -p vcf holPub.pangenie.vcf2.gz

bcftools view --threads 4 -V snps holPub.pangenie.vcf.gz \
    -r $(seq 1 29 | tr '\n' ',') \
    -Oz -o holPub.pangenie-sv.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-sv3.vcf.gz

bcftools view --threads 4 -v indels holPub.pangenie.vcf.gz \
    -r $(seq 1 29 | tr '\n' ',') \
    -Oz -o holPub.pangenie-indels.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-indels.vcf.gz

bcftools view --threads 4 -v snps holPub.pangenie.vcf.gz \
    -Oz -o holPub.pangenie-snps.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-snps.vcf.gz

####
bcftools view --threads 4 -v indels \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.indel.awk |
    bgzip -c > holPub.pangenie-indels.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-indels.vcf.gz &

bcftools view --threads 4 -v other,mnps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > holPub.pangenie-complex.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-complex.vcf.gz &

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-complex.vcf.gz \
    holPub.pangenie-indels.vcf.gz \
    -Oz -o holPub.pangenie-nonsnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-nonsnps.vcf.gz

zcat holPub.pangenie-nonsnps.vcf.gz | 
    awk '$1 ~ /#/{print; next}{split($3,a,/:|-|_/);if(a[5] > 50) {print}}' |
    bgzip -c > holPub.pangenie-sv.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-sv.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-sv.vcf.gz -Oz -o holPub.pangenie-sv.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-sv.filter.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-snps.vcf.gz -Oz -o holPub.pangenie-snps.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-snps.filter.vcf.gz

bcftools view --threads 4 -v snps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz |
    awk -f  $code/pan.snp.awk |
    bgzip -c > holPub.pangenie-snps.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-snps.vcf.gz

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-nonsnps.vcf.gz \
    holPub.pangenie-snps.vcf.gz \
    -Oz -o holPub.pangenie-var.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-var.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-var.vcf.gz -Oz -o holPub.pangenie-var.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-var.filter.vcf.gz

bcftools view --threads 4 -v mnps -r $(seq 1 29 | tr '\n' ',') holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > holPub.pangenie-mnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-mnps.vcf.gz &

l holPub.pangenie* | grep -v chr

####

zcat holPub.pangenie.filter.vcf.gz |
    awk -f $code/pangenie_remove50.awk |
    awk -f $code/pangenie.addID.awk PG="PG" |
    bgzip -c > holPub.pangenie-sv.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-sv.filter.vcf.gz

bcftools view --threads $nthreads -v snps holPub.deepv.filter.vcf.gz \
    -Oz -o holPub.deepv-snps.filter.vcf.gz && 
    tabix -f -p vcf holPub.deepv-snps.filter.vcf.gz


#####

mc_vcf=/90daydata/bull_age/liu.yang/imputation/minigraph-cactus/hol-pg2hic-2024-05-22/hol-pg2hic-2024-05-22.anno_biallelic.vcf.gz

zcat $mc_vcf |
    awk -f $code/pangenie_remove50.awk |
    awk -f $code/pangenie.addID.awk PG="PG" |
    bcftools sort -Oz -o holPub.mc-sv.vcf.gz && 
    tabix -f -p vcf holPub.mc-sv.vcf.gz

bcftools view -V snps holPub.mc-sv.vcf.gz -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    -Oz -o holPub.mc-sv3.vcf.gz && 
    tabix -f -p vcf holPub.mc-sv3.vcf.gz

bcftools view -v snps holPub.pangenie.vcf.gz \
    -Oz -o holPub.pangenie-snps.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-snps.vcf.gz

# ================================================================================================ #
# ================================================================================================ #
# ================================================================================================ #

bcftools view -h holPub.pangenie-sv.filter.vcf.gz > holPub.pangenie-sv.header
sed -i 's|#CHROM|##INFO=<ID=.,Number=1,Type=String,Description=\"aaa\">\n#CHROM|g' holPub.pangenie-sv.header

bcftools reheader -h holPub.pangenie-sv.header holPub.pangenie-sv.filter.vcf.gz -o holPub.pangenie-sv.filter2.vcf.gz

chr=10

# pangenie-sv
    bcftools view --threads $nthreads holPub.pangenie-sv.filter.vcf.gz \
        -r $chr -S holPub.samples --force-samples |
        awk -f $code/haploid2diploid.awk |
        bgzip -c -@ $nthreads > holPub.pangenie-sv.chr$chr.vcf.gz &&
        tabix -f -p vcf -@ $nthreads holPub.pangenie-sv.chr$chr.vcf.gz

# pangenie-snp
    bcftools view --threads $nthreads holPub.pangenie-snp.filter.vcf.gz \
        -r $chr -S holPub.samples --force-samples |
        awk -f $code/haploid2diploid.awk |
        bgzip -c -@ $nthreads > holPub.pangenie-snp.chr$chr.vcf.gz &&
        tabix -f -p vcf -@ $nthreads holPub.pangenie-snp.chr$chr.vcf.gz

# pangenie-var
    bcftools view --threads $nthreads holPub.pangenie-var.filter.vcf.gz \
        -r $chr -S holPub.samples --force-samples |
        awk -f $code/haploid2diploid.awk |
        bgzip -c -@ $nthreads > holPub.pangenie-var.chr$chr.vcf.gz &&
        tabix -f -p vcf -@ $nthreads holPub.pangenie-var.chr$chr.vcf.gz

## deepvariant-var
    bcftools view --threads $nthreads holPub.deepv.filter.vcf.gz \
        -r $chr -S holPub.samples --force-samples |
        awk -f $code/addID.awk PG="DV" |
        bgzip -c -@ $nthreads > holPub.deepv.chr$chr.vcf.gz

## deepvariant-snp
    bcftools view --threads $nthreads -v snps holPub.deepv.chr$chr.vcf.gz \
        -Oz -o holPub.deepv-snps.chr$chr.vcf.gz &&
        tabix -f -p vcf holPub.deepv-snps.chr$chr.vcf.gz
        
## deepvariant-indels
    bcftools view --threads $nthreads -v indels holPub.deepv.chr$chr.vcf.gz \
        -Oz -o holPub.deepv-indels.chr$chr.vcf.gz &&
        tabix -f -p vcf holPub.deepv-indels.chr$chr.vcf.gz
        
## deepvariant-mnps
    bcftools view --threads $nthreads -v mnps holPub.deepv.chr$chr.vcf.gz \
        -Oz -o holPub.deepv-mnps.chr$chr.vcf.gz &&
        tabix -f -p vcf holPub.deepv-mnps.chr$chr.vcf.gz


thin_concat holPub.pangenie-snv.chr$chr.vcf.gz "1000" "F"
thin_concat holPub.deepv.chr$chr.vcf.gz "1000" "F"
thin_concat holPub.deepv-snps.chr$chr.vcf.gz "1000 5000 10000 50000" "F"
thin_concat holPub.deepv-indels.chr$chr.vcf.gz "1000" "F"

#zcat $raw_file/$vcf_file | grep -v "##" | head -n1 | cut -f 10- | tr '\t' '\n' > run.samples.txt

# ================================================================================================ #
# ================================================================================================ #
# ================================================================================================ #
# ================================================================================================ #


proc_vcf holPub_Thin0k.pangenie-sv.chr$chr.vcf.gz
proc_vcf holPub_Thin1k.deepv-snps.chr$chr.vcf.gz 
proc_vcf holPub_Thin1k.pangenie-snv.chr$chr.vcf.gz



## thin 10k

Holstein_vcf holPub_Thin1k.deepv.chr$chr.vcf.gz

Holstein_vcf holPub_Thin1k.deepv-indels.chr$chr.vcf.gz
Holstein_vcf holPub_Thin10k.deepv-snps.chr$chr.vcf.gz 
Holstein_vcf holPub_Thin50k.deepv-snps.chr$chr.vcf.gz
Holstein_vcf holPub_Thin5k.deepv-snps.chr$chr.vcf.gz

Holstein_vcf holPub_Thin1k.pangenie-snv.chr$chr.thin10k.vcf.gz
Holstein_vcf holPub_Thin1k.pangenie-snv.chr$chr.thin5k.vcf.gz
Holstein_vcf holPub_Thin1k.pangenie-snv.chr$chr.thin50k.vcf.gz


###########
# check 

ls run_files/*gz | while read id; do 
zcat $id | 
    awk 'NR !~ "#"{
    split($3,a,/:|-|_/);
    if(a[5] > 50) {pa=1}; 
    if(a[5] <= 50) {pb=1}; 
    if(pa == pb && pa == 1){
        print "yes\t"ID;exit
    }
    }
    END{
        if(pa != pb){
            print "no\t"ID
        }
    }' ID=$id
done | awk '$1 == "no"'

for sub in pangenie pangenie-sv deepv-snps; do
   for typ in DEL INS COMPLEX; do
echo $run_file/Holstein750_Thin0k_$typ.$sub.chr$chr.vcf.gz
zcat $run_file/Holstein750_Thin0k_$typ.$sub.chr$chr.vcf.gz | awk '{split($3,a,/:|-|_/);b[a[4]]} END {for(i in b){print i}}'
done
done | less 
##########
