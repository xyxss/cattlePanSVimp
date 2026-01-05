
set -eu

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

# scp -r atlas:/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/ld_runs/holPub*/*ld* .

nt=8
proc=ld_sv.sh


# all
id=holPub.pangenie-var.filter.vcf.gz
    comp=${id/.chr10.vcf.gz}
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $PWD/$id 
"
plink --threads $nthreads --vcf ../holPub.pangenie-var.filter.vcf.gz \
--chr-set 29 --chr 1-29 --maf 0.01 \
--const-fid 0 --r2 \
--ld-window-kb 1000 \
--ld-window-r2 0.05 \
--out holPub.pangenie-var.ld \
--vcf-half-call m

awk -f $code/plink_ld_extaTypes.awk holPub.pangenie-var.ld.ld > holPub.pangenie-var.ld.ld.ext
awk '$9 !~ "snv-snv"'  holPub.pangenie-var.ld.ld.ext > holPub.pangenie-var.ld.ld.sv
awk '$1 == 10 && $9 ~ "snv-snv"'  holPub.pangenie-var.ld.ld.ext > holPub.pangenie-var.ld.ld.snp




chr=10
ids="
holPub.pangenie-sv.chr$chr.vcf.gz
holPub.pangenie-var.chr$chr.vcf.gz
holPub_Thin0k.deepv-snps.chr$chr.vcf.gz
holPub_Thin0k.pangenie-snps.chr$chr.vcf.gz
"
for id in $ids; do
    comp=${id/.chr10.vcf.gz}
    sbatch -A bull_scr \
        -D $PWD \
        --export=ALL \
        -J $comp.chr$chr.${proc//.*} \
        -c $nt \
        -o logs/$comp.chr$chr.${proc//.*}.out \
        --wrap="
bash $code/$proc $comp $chr $PWD/$id 
"
done


# all
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

# grp
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

# MAF
for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
    # MAF_vcf run_files/Holstein750_Thin0k.${sub}.chr$chr.vcf.gz
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

# typ
for sub in pangenie-var pangenie-sv deepv-snps pangenie-snps; do
    #typ_vcf run_files/Holstein750_Thin0k.${sub}.chr$chr.vcf.gz
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

# thin
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

# var1
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


# var2
for var in snps indels deepv vdeepv; do
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


