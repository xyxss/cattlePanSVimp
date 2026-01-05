slemm --lmm --phenotype_file ../data/10k.slemm.csv --bfile ../data/10k --trait QT --snp_info_file snp_info.csv --out 10k --num_threads 10
OMP_NUM_THREADS=1 slemm_gamma --pfile ../data/10k --slemm 10k --out 10k.gamma.txt
OMP_NUM_THREADS=10 slemm_gwa --pfile ../data/10k --slemm 10k --out 10k.chr1.txt --chr 1



wget https://github.com/jiang18/slemm/releases/download/v0.90.1/slemm-v0.90.1-x86_64-linux.zip

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

####
#### 5.svImp_svwgs ==================================================================================================================================


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs

nthreads=42
pre_gdata=umd50ksWGSsv_HolrePan.filter
perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv

cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld.pheno.csv
out_pre=cdcb-umd50ksv_slemm_$trait

cov=F
pc=F

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}-pc${pc}.slemm \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=8G \
        --error=${out_pre}-pc${pc}.slemm.err \
        --output=${out_pre}-pc${pc}.slemm.out \
        --time "14-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/slemm_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"

done

78886

awk 'NR == 1 || $9 < 0.05/78886' $out_pre.pc$pc.gamma.txt > $out_pre.pc$pc.gamma.sig
cat cdcb-umd50ksv_slemm_Milk.pcF.gamma.sig | 
    awk '$3 !~ "SNP" && $3 !~ "s"{print $3}' > $out_pre.pc$pc.gamma.svlist
cat $out_pre.pc$pc.gamma.sig | 
    awk '{split($3,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc /" \
    > $out_pre.pc$pc.count
    
#### 1.svImp_wgs ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs

nthreads=42
pre_gdata=umd50ksWGS_HolrePan.filter
#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm

for trait in Milk Fat Protein; do
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld.pheno.csv
out_pre=cdcb-umd50kwgsvar_slemm_$trait
cov=F
pc=F
sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}-pc${pc}.slemm \
        --cpus-per-task=$nthreads \
        --mem-per-cpu=8G \
        --error=${out_pre}-pc${pc}.slemm.err \
        --output=${out_pre}-pc${pc}.slemm.out \
        --time "14-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/slemm_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph $grm
"
done

