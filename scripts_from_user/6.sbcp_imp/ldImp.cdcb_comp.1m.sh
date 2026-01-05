
wd=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval
cd $wd
input_val=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval/All2_BTAN_WGS.pang-sv.vcf.gz
input_imp=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval/imp.pang.0.pan-snp.vcf.gz

# comp="ld-snp"

zcat $input_val | awk -f extract_SV.awk | bgzip -c > val.sv.vcf.gz && tabix -f -p vcf val.sv.vcf.gz

ln -s $input_val val.sv.vcf.gz && ln -s $input_val.tbi val.sv.vcf.gz.tbi

wd=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/1.imputation_eval
wd=/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all/2.imputation_gatk/check
nthreads=16
ls $wd/../imp.pang.*.vcf.gz |
    awk 'NR == 1' |
    while read input_imp; do
    comp=$(basename $input_imp | cut -d. -f4)
    zcat $input_imp | awk -f extract_SV.awk | bgzip -c > $comp.sv.vcf.gz && tabix -f -p vcf $comp.sv.vcf.gz
    #ln -s $input_imp $comp.sv.vcf.gz && ln -s $input_imp.tbi $comp.sv.vcf.gz.tbi
    sbatch -A bull_scr \
        -p atlas,bigmem \
        --array=0-$(($(wc -l < cv.samples.list.cv) - 1)) \
        --time=14-00:00 \
        --cpus-per-task=$nthreads \
        --job-name=2split_submit.$comp \
        --output=logs/$comp.2split_submit.%A_%a.out \
        --wrap="
bash run.split_submit.sh $comp $wd
                "
    done

cat > run.split_submit.sh << 'EOF'
set -euo pipefail

nthreads=${SLURM_CPUS_PER_TASK:-8}
samples=( $(cat cv.samples.list.cv) )
sample=${samples[$SLURM_ARRAY_TASK_ID]}
comp=$1  # pass as argument: the comparison prefix
wd=$2    # working directory

input_val=$wd/val.sv.vcf.gz

# Prepare directory
mkdir -p $wd/$comp/$sample

# Extract entire sample from imputed VCF (all chromosomes)
bcftools view --threads $nthreads -s $sample $wd/$comp.sv.vcf.gz -Oz -o $wd/$comp/$sample/imp.vcf.gz
tabix -f -p vcf $wd/$comp/$sample/imp.vcf.gz

# Extract entire sample from validation VCF (all chromosomes)
bcftools view --threads $nthreads -s $sample $input_val -Oz -o $wd/$comp/$sample/ALL.vcf.gz
tabix -f -p vcf $wd/$comp/$sample/ALL.vcf.gz

# Run downstream script
bash $wd/run.imp_cdcb.sh $comp $sample $wd 2>&1

EOF

cat > run.imp_cdcb.sh << 'EOF'
#!/bin/bash

set -euo pipefail

comp=$1
sample=$2
wd=$3
nthreads=$SLURM_CPUS_PER_TASK

ref_path=/90daydata/bull_age/liu.yang/ref
ref_fa=$ref_path/ARS_UCD_v2.0.fa
ref_rm=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat
ref_gff=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_gff

function main () {
    cd $wd/$comp/$sample
    bcftools view --threads $nthreads -i 'GT="0/0"' ALL.vcf.gz -Oz -o REF.vcf.gz &&
    tabix -f -p vcf REF.vcf.gz
    bcftools view --threads $nthreads -i 'GT="1/1" || GT~"0/1" || GT~"1/0"' ALL.vcf.gz -Oz -o ALT.vcf.gz &&
    tabix -f -p vcf ALT.vcf.gz
    bcftools view --threads $nthreads -i 'GT="1/1"' ALL.vcf.gz -Oz -o HOM.vcf.gz &&
    tabix -f -p vcf HOM.vcf.gz
    bcftools view --threads $nthreads -i 'GT~"0/1" || GT~"1/0"' ALL.vcf.gz -Oz -o HET.vcf.gz &&
    tabix -f -p vcf HET.vcf.gz

    gts="ALL REF ALT HET HOM"
    for gt in $gts; do
        g_vailFile_forGT $gt
    done

    rms="None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"
    for rm in $rms; do
        g_vailFile_forRM $ref_rm/rm.$rm.bed $rm
    done

    gfs="nongene intron antisense_RNA CDS C_gene_segment D_loop exon gene lnc_RNA miRNA mRNA ncRNA primary_transcript pseudogene RNase_MRP_RNA rRNA snoRNA snRNA SRP_RNA telomerase_RNA transcript tRNA V_gene_segment"
    for gf in $gfs; do
        g_vailFile_forRM $ref_gff/gff.$gf.bed $gf
    done

    rm -rf $wd/$comp/$sample
}

function g_vailFile_forGT () {
    local $gt=$1

    bcftools query -f '%CHROM\t%POS\n' $gt.vcf.gz > $gt.snp_list.txt
    bcftools view --threads $nthreads -R $gt.snp_list.txt imp.vcf.gz -Oz -o $gt.imp.in.vcf.gz &&
    tabix -f -p vcf $gt.imp.in.vcf.gz

    rtg vcfeval -b $gt.vcf.gz \
        -c $gt.imp.in.vcf.gz \
        --Xmax-length=1000000 \
        -t ${ref_fa/fa/sdf} -T $nthreads --no-roc \
        -o $gt
    
    cat $gt/summary.txt | 
        awk '$1 == "None"{printf ID" ";for(i=1;i<=NF;i++){printf $i" "};printf "\n"}' \
            ID="$comp $sample $gt" >> $wd/$comp.table2.txt &&
    rm -rf $gt
}

function g_vailFile_forRM () {
    local bed=$1
    local rm=$2

    rtg vcfeval -b ALL.vcf.gz -c imp.vcf.gz \
        --Xmax-length=1000000 \
        -e $bed -t ${ref_fa/fa/sdf} -T $nthreads --no-roc -o $rm

    cat $rm/summary.txt | 
        awk '$1 == "None"{printf ID" ";for(i=1;i<=NF;i++){printf $i" "};printf "\n"}' \
            ID="$comp $sample $rm" >> $wd/$comp.table2.txt &&
    rm -rf $rm
}

main
EOF



### 


/90daydata/bull_age/liu.yang/imputation/pangenie_cdcb/3.pan_all

ls */check/*gz | while read id; do
    iid=$(awk -F "[_/.]" '{print $3"."$5}' <<< $id)
    zcat $id | awk ' $1 !~ "^#"{split($8,a, /[;=]/);print $3,a[2],a[4],a[7],a[9]}' > $iid.bed
    done

C:\Users\Liu.Yang\OneDrive - University of Maryland\Data\2024-02-07.cattlePanSV\9.imp_panel
