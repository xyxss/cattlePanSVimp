cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/slemm_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6

nthreads=$SLURM_CPUS_PER_TASK
slemm_path=/project/bull_scr/liu.yang/software/slemm-v0.90.1-x86_64-linux/

if [[ $cov == "F" ]]; then
    $slemm_path/slemm --num_threads $nthreads --bfile $pre_gdata \
        --max_heritability 0.4 \
        --phenotype_file $ph --trait $trait \
        --lmm --snp_info_file $pre_gdata.snp_info.csv \
        --min_maf 0.01 \
        --seed 9 \
        --output_file $out_pre.pc$pc
else
    $slemm_path/slemm --num_threads $nthreads --bfile $pre_gdata \
        --max_heritability 0.4 \
        --phenotype_file $ph --trait $trait \
        --lmm --snp_info_file $pre_gdata.snp_info.csv \
        --covariate_file $cov \
        --covariate_names $(printf "PC%d," $(seq 1 $pc)) \
        --min_maf 0.01 \
        --seed 9 \
        --output_file $out_pre.pc$pc
fi

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \
    $slemm_path/slemm_gamma \
    --pfile $pre_gdata \
    --slemm $out_pre.pc$pc \
    --out $out_pre.pc$pc.gamma.txt

for chr in $(seq 1 29); do
    sbatch -A bull_scr --partition atlas \
            --job-name=${out_pre}-pc${pc}.slemm_gwa.chr$chr \
            --cpus-per-task=$nthreads \
            --error=${out_pre}-pc${pc}.slemm_gwa.chr$chr.err \
            --output=${out_pre}-pc${pc}.slemm_gwa.chr$chr.out \
            --time "14-00:00:00" \
            --wrap="
OMP_NUM_THREADS=$nthreads \
    $slemm_path/slemm_gwa \
    --pfile $pre_gdata \
    --slemm $out_pre.pc$pc \
    -c $chr \
    --out $out_pre.pc$pc.gwa.chr$chr.txt
        "
done
EOF


#### 5.svImp_svwgs ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85

plink2 --threads 8 --bfile umd50ksWGS_HolrePan.filter.r85 --chr-set 29 \
    --make-pgen --out umd50ksWGS_HolrePan.filter.r85


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85

nthreads=42
pre_gdata=umd50ksWGS_HolrePan.filter.r85
#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cat $pre_gdata.pvar | awk 'BEGIN{print "SNP"} $0 !~ "#"{print $3}' > $pre_gdata.snp_info.csv

grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm
cov=F
pc=F
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}-pc${pc}.slemm \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=8G \
        --error=logs/${out_pre}-pc${pc}.slemm.err \
        --output=logs/${out_pre}-pc${pc}.slemm.out \
        --time "14-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/slemm_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

done

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait

awk 'NR == 1 || $9 < 0.05/6048794' $out_pre.pc$pc.gwa.chr*txt \
    >  $out_pre.pc$pc.gwa.sig
cat $out_pre.pc$pc.gwa.sig | 
    awk '$3 !~ "SNP" && $3 !~ "s"{print}' > $out_pre.pc$pc.svlist
cat $out_pre.pc$pc.gwa.sig | 
    awk '{split($3,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc /" \
    > $out_pre.pc$pc.svlist.count

done


### qq

thr=0.05/6048794
thr=8.26611e-09
for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk '($9 < Thr*1e1 && NR % 2 == 1) || ($9 < Thr*1e2 && NR % 5 == 1) || ($9 < Thr*1e3 && NR % 10 == 1) || ($9 > Thr*1e3 && NR % 200 == 1)' Thr=$thr > sigplot.$out_pre.pc$pc.gwa.txt

cat $out_pre.pc$pc.gamma.txt | awk '($9 < Thr*1e1 && NR % 2 == 1) || ($9 < Thr*1e2 && NR % 5 == 1) || ($9 < Thr*1e3 && NR % 10 == 1) || ($9 > Thr*1e3 && NR % 200 == 1)' Thr=$thr > sigplot.$out_pre.pc$pc.gamma.txt
done

scp atlas:/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85/*gwa.sig .
scp atlas:/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85/sigplot*gwa.txt .
scp atlas:/90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/0.svImp_r85/*.svlist.count .



cdcb-umd50ksv_slemm_Body_depth.pcF  1
cdcb-umd50ksv_slemm_Body_depth.pcF SNP 664
cdcb-umd50ksv_slemm_Body_depth.pcF DEL 3
cdcb-umd50ksv_slemm_Body_depth.pcF INS 3
cdcb-umd50ksv_slemm_Body_depth.pcF sDEL 153
cdcb-umd50ksv_slemm_Body_depth.pcF sINS 173
cdcb-umd50ksv_slemm_Body_depth.pcF COMPLEX 4
cdcb-umd50ksv_slemm_Body_depth.pcF sCOMPLEX 67
cdcb-umd50ksv_slemm_Fat.pcF  1
cdcb-umd50ksv_slemm_Fat.pcF SNP 2241
cdcb-umd50ksv_slemm_Fat.pcF DEL 10
cdcb-umd50ksv_slemm_Fat.pcF INS 19
cdcb-umd50ksv_slemm_Fat.pcF sDEL 465
cdcb-umd50ksv_slemm_Fat.pcF sINS 445
cdcb-umd50ksv_slemm_Fat.pcF COMPLEX 9
cdcb-umd50ksv_slemm_Fat.pcF sCOMPLEX 218
cdcb-umd50ksv_slemm_Milk.pcF  1
cdcb-umd50ksv_slemm_Milk.pcF SNP 1494
cdcb-umd50ksv_slemm_Milk.pcF DEL 6
cdcb-umd50ksv_slemm_Milk.pcF INS 17
cdcb-umd50ksv_slemm_Milk.pcF sDEL 342
cdcb-umd50ksv_slemm_Milk.pcF sINS 373
cdcb-umd50ksv_slemm_Milk.pcF COMPLEX 11
cdcb-umd50ksv_slemm_Milk.pcF sCOMPLEX 147
cdcb-umd50ksv_slemm_Protein.pcF  1
cdcb-umd50ksv_slemm_Protein.pcF SNP 825
cdcb-umd50ksv_slemm_Protein.pcF DEL 2
cdcb-umd50ksv_slemm_Protein.pcF INS 7
cdcb-umd50ksv_slemm_Protein.pcF sDEL 153
cdcb-umd50ksv_slemm_Protein.pcF sINS 175
cdcb-umd50ksv_slemm_Protein.pcF COMPLEX 2
cdcb-umd50ksv_slemm_Protein.pcF sCOMPLEX 63
cdcb-umd50ksv_slemm_Stature.pcF  1
cdcb-umd50ksv_slemm_Stature.pcF SNP 978
cdcb-umd50ksv_slemm_Stature.pcF DEL 6
cdcb-umd50ksv_slemm_Stature.pcF INS 7
cdcb-umd50ksv_slemm_Stature.pcF sDEL 264
cdcb-umd50ksv_slemm_Stature.pcF sINS 293
cdcb-umd50ksv_slemm_Stature.pcF COMPLEX 3
cdcb-umd50ksv_slemm_Stature.pcF sCOMPLEX 132


thr=8.26611e-09

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50ksv_slemm_$trait
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt
done

