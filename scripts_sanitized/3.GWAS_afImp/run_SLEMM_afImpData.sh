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
    sbatch \
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


#### Imp GWAS for SV only ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/0.svImp_r8

awk '$3 !~ "SNP" && $3 !~ "s" && $7 ~ "IMP" {split($7,a,/;|=/);if(a[3] > 0.8)print} ' \
    ../genotypes/umdsample.svwgssnp.afterImput.filter.pvar > umd50ksWGS_HolrePan.filter.r8.sv

plink2 --threads 8 --pfile ../genotypes/umdsample.svwgssnp.afterImput.filter --chr-set 29 \
    --extract umd50ksWGS_HolrePan.filter.r8.sv --make-pgen --out umd50ksWGS_HolrePan.filter.r8sv

plink2 --threads 8 --pfile ../genotypes/umdsample.svwgssnp.afterImput.filter --chr-set 29 \
    --extract umd50ksWGS_HolrePan.filter.r8.sv --make-bed --out umd50ksWGS_HolrePan.filter.r8sv


nthreads=42
pre_gdata=umd50ksWGS_HolrePan.filter.r8sv
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
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt
done


#### Imp GWAS for SV&SNP ==================================================================================================================================

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/0.svImp_r8

awk '$3 !~ "SNP" && $3 !~ "s" && $7 ~ "IMP" {split($7,a,/;|=/);if(a[3] > 0.8)print} ' \
    ../genotypes/umdsample.svwgssnp.afterImput.filter.pvar > umd50ksWGS_HolrePan.filter.r8.sv

awk '{if($7 !~ "IMP"){print} else {split($7,a,/;|=/);if(a[3] > 0.8)print}} ' \
    ../genotypes/umdsample.svwgssnp.afterImput.filter.pvar > umd50ksWGS_HolrePan.filter.r8

plink2 --threads 8 --pfile ../genotypes/umdsample.svwgssnp.afterImput.filter --chr-set 29 \
    --extract umd50ksWGS_HolrePan.filter.r8 --make-pgen --out umd50ksWGS_HolrePan.filter.r8

plink2 --threads 8 --pfile ../genotypes/umdsample.svwgssnp.afterImput.filter --chr-set 29 \
    --extract umd50ksWGS_HolrePan.filter.r8 --make-bed --out umd50ksWGS_HolrePan.filter.r8


nthreads=42
pre_gdata=umd50ksWGS_HolrePan.filter.r8
#perl -e 'print "SNP\n"; while(<>) {@c=split /\s+/; print "$c[1]\n"}' < $pre_gdata.bim > $pre_gdata.snp_info.csv
cat $pre_gdata.pvar | awk 'BEGIN{print "SNP"} $0 !~ "#"{print $3}' > $pre_gdata.snp_info.csv

grm=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.chip.grm
cov=F
pc=F
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/hol.yld_type.pheno.csv

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50k_slemm_$trait

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

#
for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50k_slemm_$trait
pc=F
cat $out_pre.pc$pc.gwa.chr*txt | awk 'NR == 1 || $1 !~ "#"' > $out_pre.pc$pc.gwa.txt && 
rm $out_pre.pc$pc.gwa.chr*txt
done

for trait in Milk Fat Protein Stature Body_depth; do
out_pre=cdcb-umd50k_slemm_$trait

awk 'NR == 1 || $11 < 0.05/10713446' $out_pre.pc$pc.gwa.txt \
    >  $out_pre.pc$pc.gwa.sig
cat $out_pre.pc$pc.gwa.sig | 
    awk '$3 !~ "SNP" && $3 !~ "s"{print}' > $out_pre.pc$pc.svlist
cat $out_pre.pc$pc.gwa.sig | 
    awk '{split($3,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc /" \
    > $out_pre.pc$pc.svlist.count
done

