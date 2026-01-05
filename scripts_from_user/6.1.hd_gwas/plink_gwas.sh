
#cat > /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh <<"EOF"
pre_gdata=$1
out_pre=$2
trait=$3
cov=$4
pc=$5
ph=$6

chrom=29
nthreads=$SLURM_CPUS_PER_TASK

if [[ $cov == "F" ]]; then
    plink2 --threads $nthreads --bfile $pre_gdata \
        --chr-set $chrom \
        --pheno iid-only $ph --pheno-name $trait \
        --no-categorical --glm allow-no-covars \
        --no-input-missing-phenotype \
        --out $out_pre.pc$pc
else
    plink2 --threads $nthreads --bfile $pre_gdata --chr-set 29 \
        --pheno iid-only $ph --pheno-name $trait \
        --covar $cov --covar-name $(printf "PC%d " $(seq 1 $pc)) \
        --no-categorical --glm hide-covar \
        --no-input-missing-phenotype \
        --out $out_pre.pc$pc
fi
flt=$(awk '$0 ~ "variants loaded from"{print $1}' $out_pre.pc$pc.log)
awk '$15 < 0.05/flt && $10 == "ADD"{print $3, $7, $4, $9, $12, $13, $15, $11}' \
    flt=$flt $out_pre.pc$pc.$trait.glm.linear > $out_pre.pc$pc.$trait.summary_stats.txt
cat $out_pre.pc$pc.$trait.summary_stats.txt | 
    awk '$1 !~ "SNP" && $1 !~ "s"{print $1}' > $out_pre.pc$pc.$trait.svlist
cat $out_pre.pc$pc.$trait.summary_stats.txt | 
    awk '{split($1,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}' |
    sed "s/^/$out_pre.pc$pc.$trait /" \
    > $out_pre.pc$pc.$trait.svlist.count
EOF


nthreads=16
for var in sv var gsnp; do

pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
out_pre=cdcb_gwas_nocv.$var
trait=Milk
cov=All2_BTAN_WGS.var.bfile.maf1.eigenvec
pc="2"
ph=ph.hol.traits

for trait in fat milk pro; do

sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done
done

nthreads=16
for var in sv var gsnp; do

pre_gdata=All2_BTAN_WGS.$var.bfile.maf1
out_pre=cdcb_gwas_nocv.$var
trait=Milk
cov=F
pc=F
ph=cdcb.hol.traits

for trait in fat milk pro; do

sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done
done



cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/1.svImp_wgs

nthreads=48
pre_gdata=umd50ksWGS_HolrePan.filter
out_pre=cdcb-umd50k_gwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=2
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done

cat *.svlist.count
cdcb-umd50k_gwas.pc10.Fat SNP 3404931
cdcb-umd50k_gwas.pc10.Fat DEL 5001
cdcb-umd50k_gwas.pc10.Fat INS 8660
cdcb-umd50k_gwas.pc10.Fat sDEL 290739
cdcb-umd50k_gwas.pc10.Fat sINS 328239
cdcb-umd50k_gwas.pc10.Fat COMPLEX 7925
cdcb-umd50k_gwas.pc10.Fat sCOMPLEX 129062
cdcb-umd50k_gwas.pc10.Milk SNP 4274987
cdcb-umd50k_gwas.pc10.Milk DEL 6265
cdcb-umd50k_gwas.pc10.Milk sDEL 364881
cdcb-umd50k_gwas.pc10.Milk INS 10783
cdcb-umd50k_gwas.pc10.Milk sINS 412991
cdcb-umd50k_gwas.pc10.Milk COMPLEX 9708
cdcb-umd50k_gwas.pc10.Milk sCOMPLEX 165298
cdcb-umd50k_gwas.pc10.Protein SNP 4825174
cdcb-umd50k_gwas.pc10.Protein DEL 7044
cdcb-umd50k_gwas.pc10.Protein sDEL 411631
cdcb-umd50k_gwas.pc10.Protein INS 12212
cdcb-umd50k_gwas.pc10.Protein sINS 467500
cdcb-umd50k_gwas.pc10.Protein COMPLEX 11157
cdcb-umd50k_gwas.pc10.Protein sCOMPLEX 182702
cdcb-umd50k_gwas.pc2.Fat SNP 4923590
cdcb-umd50k_gwas.pc2.Fat DEL 7271
cdcb-umd50k_gwas.pc2.Fat sDEL 421288
cdcb-umd50k_gwas.pc2.Fat INS 12677
cdcb-umd50k_gwas.pc2.Fat sINS 477996
cdcb-umd50k_gwas.pc2.Fat COMPLEX 11535
cdcb-umd50k_gwas.pc2.Fat sCOMPLEX 188906
cdcb-umd50k_gwas.pc2.Milk SNP 5474491
cdcb-umd50k_gwas.pc2.Milk DEL 8070
cdcb-umd50k_gwas.pc2.Milk INS 13926
cdcb-umd50k_gwas.pc2.Milk sDEL 467468
cdcb-umd50k_gwas.pc2.Milk sINS 530885
cdcb-umd50k_gwas.pc2.Milk COMPLEX 12739
cdcb-umd50k_gwas.pc2.Milk sCOMPLEX 209245
cdcb-umd50k_gwas.pc2.Protein SNP 5845628
cdcb-umd50k_gwas.pc2.Protein DEL 8553
cdcb-umd50k_gwas.pc2.Protein sDEL 499254
cdcb-umd50k_gwas.pc2.Protein INS 14926
cdcb-umd50k_gwas.pc2.Protein sINS 566983
cdcb-umd50k_gwas.pc2.Protein COMPLEX 13558
cdcb-umd50k_gwas.pc2.Protein sCOMPLEX 223530


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/4.snp.gwas
nthreads=48
pre_gdata=../allseq_1kbulls.hol.chrall.maf01
out_pre=cdcb-umd50k_snpgwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Pro_Percent Milk Fat_Percent Fat Protein; do
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/2.svImp_ld

nthreads=48
pre_gdata=umd50ksLD_HolrePan.filter2
out_pre=cdcb-umd50k_ldgwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=2
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done



cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

done
done


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do
sbatch -A bull_scr \
        --job-name=${out_pre}_${trait}.GWAS \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GWAS.err \
        --output=${out_pre}_${trait}.GWAS.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

done
done


cd /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/5.svImp_svwgs

nthreads=48
pre_gdata=umd50ksWGSsv_HolrePan.filter
out_pre=cdcb-umd50ksv_wgsgwas
#trait=Milk
cov=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.chip.maf01.eigenvec
pc=10
ph=/90daydata/bull_age/liu.yang/imputation/hd_gwas/pheotypes/hol.yld.pheno.phe2

for trait in Milk Fat Protein; do

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"

cov=F
pc=F

sbatch -A bull_scr --partition atlas \
        --job-name=${out_pre}_${trait}.GCTA \
        --cpus-per-task=$nthreads \
        --error=${out_pre}_${trait}.GCTA.err \
        --output=${out_pre}_${trait}.GCTA.out \
        --time "1-00:00:00" \
        --wrap="
bash /90daydata/bull_age/liu.yang/imputation/hd_gwas/genotypes/plink_gwas.sh \
    $pre_gdata $out_pre $trait $cov $pc $ph
"
done

cat *.svlist.count
cdcb-umd50ksv_wgsgwas.pc10.Fat DEL 6506
cdcb-umd50ksv_wgsgwas.pc10.Fat INS 11242
cdcb-umd50ksv_wgsgwas.pc10.Fat COMPLEX 10348
cdcb-umd50ksv_wgsgwas.pc10.Milk DEL 7694
cdcb-umd50ksv_wgsgwas.pc10.Milk INS 13342
cdcb-umd50ksv_wgsgwas.pc10.Milk COMPLEX 12016
cdcb-umd50ksv_wgsgwas.pc10.Protein DEL 8498
cdcb-umd50ksv_wgsgwas.pc10.Protein INS 14691
cdcb-umd50ksv_wgsgwas.pc10.Protein COMPLEX 13497
cdcb-umd50ksv_wgsgwas.pcF.Fat DEL 14886
cdcb-umd50ksv_wgsgwas.pcF.Fat INS 26301
cdcb-umd50ksv_wgsgwas.pcF.Fat COMPLEX 23924
cdcb-umd50ksv_wgsgwas.pcF.Milk DEL 14544
cdcb-umd50ksv_wgsgwas.pcF.Milk INS 25641
cdcb-umd50ksv_wgsgwas.pcF.Milk COMPLEX 23309
cdcb-umd50ksv_wgsgwas.pcF.Protein DEL 15028
cdcb-umd50ksv_wgsgwas.pcF.Protein INS 26482
cdcb-umd50ksv_wgsgwas.pcF.Protein COMPLEX 24177