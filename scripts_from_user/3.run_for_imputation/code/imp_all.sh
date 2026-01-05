#!/bin/bash

#set -o nounset
#set -o errexit

comp=$1
chr=$2
input_vcf_gz=$3
options=$4
nthreads=$SLURM_CPUS_PER_TASK
#pangenieSV
rtgNcpu=1


####
work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/imp_runs
code=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/00.codes

function main (){

    options="${options#options=}"
    # Split the options by comma
    IFS=',' read -r -a option_pairs <<< "$options"

    declare -A optionDict
    for pair in "${option_pairs[@]}"; do
        IFS=':' read -r key value <<< "$pair"
        optionDict["$key"]="$value"
    done

## 
    ref_path=/90daydata/bull_age/liu.yang/ref
    ref_fa=$ref_path/ARS_UCD_v2.0.fa
    ref_rm=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat
    ref_gff=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_gff

    minimac4=/home/liu.yang/software/minimac4-4.1.6-Linux-x86_64/bin/minimac4
    beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar
    #eagle=/home/liu.yang/software/Eagle_v2.4.1/eagle

    source /home/liu.yang/miniconda/etc/profile.d/conda.sh
    
    mkdir -p $work_dir/$comp/$comp.chr$chr
    cd $work_dir/$comp/$comp.chr$chr

    f="cv.samples.list.cv"
    if [[ ! -f $f || $(stat -c%s "$f") -eq 0 ]]; then
        zcat $input_vcf_gz | grep -v "##" | head -n 1 | 
            awk -f $code/cv_samples.awk > cv.samples.list.cv
    fi

    split_train_vali_samples

    ref2Beagle5
    if [ "${optionDict["tool"]}" == "T" ]; then
        ref2minimac4
    fi

    if [ "${optionDict["tool"]}" == "T" ]; then
        ref2minimac4
    fi

    if [ "${optionDict["rate"]}" != "F" ]; then
        rates="0.1 0.2 0.5"
        for rate in $rates ;do
            val_misrate
        done

        if [ "${optionDict["rate"]}" == "T" ]; then
            rate="0"
            val_misvar
            rates="0 0.1 0.2 0.5"
        fi

        if [ "${optionDict["rate"]}" == "X" ]; then
            rates="0.1 0.2 0.5"
        fi

    else
        rates="0"
        val_misvar
    fi

    for rate in $rates ;do
        imp2Beagle5
        if [ "${optionDict["tool"]}" == "T" ]; then
            imp2minimac4
        fi
    done
    
    
    echo "Comp Rate Chr Sample GT Tool Threshold True-pos-baseline True-pos-call False-pos False-neg Precision Sensitivity F-measure" > $work_dir/$comp/$comp.chr$chr.check.table2.GT.txt
    if [ "${optionDict["rm"]}" == "T" ]; then
        echo "Comp Rate Chr Sample GT Tool Threshold True-pos-baseline True-pos-call False-pos False-neg Precision Sensitivity F-measure" > $work_dir/$comp/$comp.chr$chr.check.table2.RM.txt
    fi

    for rate in $rates ;do
        mkdir -p set.missing$rate

        cat cv.samples.list.cv | while read sample ; do
            g_vailFile
        done
    done

}

function sample_val () {
    mkdir -p set.missing$rate/$sample/set.IDs.$gt.dir

    awk 'FNR==NR{a[$1]++;next} (a[$3] || $1 ~ "#"){print}' \
        set.missing$rate/$sample/set.IDs.$gt <(bcftools view -s $sample validation.vcf.gz) |
        bgzip -c > set.missing$rate/$sample/set.IDs.$gt.dir/validation.vcf.gz && 
        tabix -p vcf set.missing$rate/$sample/set.IDs.$gt.dir/validation.vcf.gz
}

function g_vailFile_forGT () {
    local tool=$1
    mkdir -p set.missing$rate/$sample/set.IDs.$gt.dir
    awk 'FNR==NR{a[$1]++;next} (a[$3] || $1 ~ "#"){print}' \
        set.missing$rate/$sample/set.IDs.$gt \
        <(bcftools view -s $sample set.missing$rate/validation.$tool.miss.vcf.gz) |
        bgzip -c > set.missing$rate/$sample/set.IDs.$gt.dir/$tool.vcf.gz && 
        tabix -p vcf set.missing$rate/$sample/set.IDs.$gt.dir/$tool.vcf.gz

    rtg vcfeval -b set.missing$rate/$sample/set.IDs.$gt.dir/validation.vcf.gz \
        -c set.missing$rate/$sample/set.IDs.$gt.dir/$tool.vcf.gz \
        -t ${ref_fa/fa/sdf} -T $rtgNcpu --no-roc --Xmax-length=1000000 \
        -o set.missing$rate/$sample/set.IDs.$gt.dir/$tool.$sample.mIDs
    
    cat set.missing$rate/$sample/set.IDs.$gt.dir/$tool.$sample.mIDs/summary.txt | 
        awk '$1 == "None"{printf ID" ";for(i=1;i<=NF;i++){printf $i" "};printf "\n"}' \
            ID="$comp $rate $chr $sample $gt $tool" >> $work_dir/$comp/$comp.chr$chr.check.table2.GT.txt &&
    rm -rf set.missing$rate/$sample/set.IDs.$gt.dir/$tool.$sample.mIDs
}

function g_vailFile_forRM () {
    local tool=$1
    local bed=$2
    local rm=$3
    rtg vcfeval -b set.missing$rate/$sample/set.IDs.ALT.dir/validation.vcf.gz \
        -c set.missing$rate/$sample/set.IDs.ALT.dir/$tool.vcf.gz \
        -e $bed -t ${ref_fa/fa/sdf} -T $rtgNcpu --no-roc --Xmax-length=1000000 \
        -o set.missing$rate/$sample/set.IDs.ALT.dir/$tool.$sample.mIDs.$rm
    
    cat set.missing$rate/$sample/set.IDs.ALT.dir/$tool.$sample.mIDs.$rm/summary.txt | 
        awk '$1 == "None"{printf ID" ";for(i=1;i<=NF;i++){printf $i" "};printf "\n"}' \
            ID="$comp $rate $chr $sample $rm $tool" >> $work_dir/$comp/$comp.chr$chr.check.table2.RM.txt &&
    rm -rf set.missing$rate/$sample/set.IDs.ALT.dir/$tool.$sample.mIDs.$rm
}

function g_vailFile () {

    count=0

    mkdir -p set.missing$rate/$sample/
    cp set.missing$rate/set.$sample.missing$rate.IDs set.missing$rate/$sample/set.IDs.ALL
    awk '$2="0/0"' set.missing$rate/set.$sample.missing$rate.IDs > set.missing$rate/$sample/set.IDs.REF
    awk '$2!="0/0"' set.missing$rate/set.$sample.missing$rate.IDs > set.missing$rate/$sample/set.IDs.ALT
    awk '$2=="1/1"' set.missing$rate/set.$sample.missing$rate.IDs > set.missing$rate/$sample/set.IDs.HOM
    awk '$2=="0/1" || $2=="1/0"' set.missing$rate/set.$sample.missing$rate.IDs > set.missing$rate/$sample/set.IDs.HET

    if [ "${optionDict["gt"]}" == "T" ]; then
        gts="ALL REF ALT HET HOM"
    else
        gts="ALT"
    fi
    for gt in $gts; do
        sample_val

        if [ "${optionDict["tool"]}" == "T" ]; then
            g_vailFile_forGT Beagle5 && g_vailFile_forGT Minimac4 &
        else
            g_vailFile_forGT Beagle5 &
        fi

        ((count++))
        if (( count % nthreads == 0 )); then
            wait
        fi
    done
    

    if [ "${optionDict["typ"]}" == "T" ]; then
        types="$(awk '{split($1,a,/:|-|_/);b[a[4]]} END {for(i in b){print i}}' set.missing$rate/$sample/set.IDs.ALT)"
        for typ in $types; do
            awk '$1~typ' typ=$typ set.missing$rate/$sample/set.IDs.ALT > set.missing$rate/$sample/set.IDs.$typ
            gt=$typ
            sample_val
            if [ "${optionDict["tool"]}" == "T" ]; then
                g_vailFile_forGT Beagle5 && g_vailFile_forGT Minimac4 &
            else
                g_vailFile_forGT Beagle5 &
            fi

            ((count++))
            if (( count % nthreads == 0 )); then
                wait
            fi
        done
    fi

    if [ "${optionDict["rm"]}" == "T" ]; then
        rms="None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"
        for rm in $rms; do
            g_vailFile_forRM Beagle5 $ref_rm/rm.$rm.bed $rm && 
            g_vailFile_forRM Minimac4 $ref_rm/rm.$rm.bed $rm &

            ((count++))
            if (( count % nthreads == 0 )); then
                wait
            fi
        done
    fi

    if [ "${optionDict["gf"]}" == "T" ]; then
        gfs="nongene
intron
antisense_RNA
CDS
C_gene_segment
D_loop
exon
gene
lnc_RNA
miRNA
mRNA
ncRNA
primary_transcript
pseudogene
RNase_MRP_RNA
rRNA
snoRNA
snRNA
SRP_RNA
telomerase_RNA
transcript
tRNA
V_gene_segment"
        for gf in $gfs; do
            g_vailFile_forRM Beagle5 $ref_gff/gff.$gf.bed $gf && 
            g_vailFile_forRM Minimac4 $ref_gff/gff.$gf.bed $gf &

            ((count++))
            if (( count % nthreads == 0 )); then
                wait
            fi
        done
    fi



    wait
}

function ref2Beagle5() {
    java -jar $beagle gt=training.vcf.gz out=reference
}

function imp2Beagle5() {
    java -jar $beagle ref=reference.vcf.gz \
        gt=set.missing$rate/validation.missing.vcf.gz \
        out=set.missing$rate/validation.Beagle.miss

    zcat set.missing$rate/validation.Beagle.miss.vcf.gz | 
        awk '{if($0 ~ "##source"){$0=$0"\n##contig=<ID="chr">"}  print}' chr=$chr | 
        bgzip -c > set.missing$rate/validation.Beagle5.miss.vcf.gz
}

function ref2minimac4() {
    zcat reference.vcf.gz | 
        awk '{if($0 ~ "##source"){$0=$0"\n##contig=<ID="chr">"} print}' chr=$chr | 
        bgzip -c > reference1.vcf.gz

    $minimac4 --threads $nthreads \
        --compress-reference reference1.vcf.gz \
        > reference.msav
}

function imp2minimac4() {
    $minimac4 --threads $nthreads \
        reference.msav \
        set.missing$rate/validation.missing.vcf.gz \
        --format GT \
        -O vcf.gz \
        -o set.missing$rate/validation.Minimac4.miss.vcf.gz
}

function val_misrate () {
    mkdir -p set.missing$rate

    bcftools view --threads $nthreads validation.vcf.gz | awk -f $code/set_vcf_ID.missing.awk rate=$rate |
    bgzip -c > set.missing$rate/validation.missing.vcf.gz && \
    tabix -p vcf set.missing$rate/validation.missing.vcf.gz
}

function val_misvar () {
    mkdir -p set.missing$rate

    bcftools view --threads $nthreads validation.vcf.gz | awk -f $code/set_vcf_SV.missing.awk rate=$rate |
        bgzip -c > set.missing$rate/validation.missing.vcf.gz && \
        tabix -p vcf set.missing$rate/validation.missing.vcf.gz
}

function split_train_vali_samples () {
    bcftools view --threads $nthreads -S cv.samples.list.cv $input_vcf_gz --force-samples | bgzip -c > validation.vcf.gz
    bcftools view --threads $nthreads -S ^cv.samples.list.cv $input_vcf_gz --force-samples | bgzip -c > training.vcf.gz
}

main


