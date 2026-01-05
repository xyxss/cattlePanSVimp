
ref_path=/90daydata/bull_age/liu.yang/ref
ref_fa=$ref_path/ARS_UCD_v2.0.fa
ref_rm=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat


#### data preparation
work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp
mkdir -p $work_dir/pca_runs && cd $work_dir/pca_runs

run_file=$work_dir/run_files
code=$work_dir/00.codes


nthreads=$SLURM_CPUS_PER_TASK

### pangenie-sv run plink pca

plink2 --vcf ../holPub.pangenie-sv.vcf.gz \
    --chr-set 29 --chr 1-29 --maf 0.01 --pca \
    --out holPub.pangenie-sv.pca_results \
    --threads $nthreads

plink2 --threads $nthreads \
    --vcf ../holPub.pangenie-var.filter.vcf.gz \
    --pca --chr-set 29 --chr 1-29 --maf 0.01 \
    --vcf-half-call m \
    --out holPub.pangenie-var.pca_results

plink2 --threads $nthreads \
    --vcf ../holPub.pangenie-snps.filter.vcf.gz \
    --pca --chr-set 29 --chr 1-29 --maf 0.01 \
    --vcf-half-call m \
    --out holPub.pangenie-snps.pca_results
    
plink2 --threads $nthreads \
    --vcf ../holPub.deepv.filter.vcf.gz \
    --pca --chr-set 29 --chr 1-29 --maf 0.01 \
    --vcf-half-call m \
    --out holPub.deepv.pca_results

plink2 --threads $nthreads \
    --vcf ../holPub.deepv-snps.filter.vcf.gz \
    --pca --chr-set 29 --chr 1-29 --maf 0.01 \
    --vcf-half-call m \
    --out holPub.deepv.pca_results

cat ../holPub.samples.group | awk -F "\t" '{print $1 > $2}' 

vcftools --gzvcf ../holPub.pangenie-sv.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Jersey \
    --maf 0.01 \
    --out Holstein-Jersey.fst &

vcftools --gzvcf ../holPub.pangenie-sv.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Holstein-X-Jersey \
    --maf 0.01 \
    --out Holstein-Cross.fst &

vcftools --gzvcf ../holPub.pangenie-sv.filter.vcf.gz \
    --weir-fst-pop Holstein-X-Jersey \
    --weir-fst-pop Jersey \
    --maf 0.01 \
    --out Cross-Jersey.fst &

# pangenie-snps
vcftools --gzvcf ../holPub.pangenie-snps.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Holstein-Jersey.pangenie-snps.fst &

vcftools --gzvcf ../holPub.pangenie-snps.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Holstein-X-Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Holstein-Cross.pangenie-snps.fst &

vcftools --gzvcf ../holPub.pangenie-snps.filter.vcf.gz \
    --weir-fst-pop Holstein-X-Jersey \
    --weir-fst-pop Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Cross-Jersey.pangenie-snps.fst &

# deepv-snp
vcftools --gzvcf ../holPub.deepv-snps.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Holstein-Jersey.deepv-snps.fst &

vcftools --gzvcf ../holPub.deepv-snps.filter.vcf.gz \
    --weir-fst-pop Holstein \
    --weir-fst-pop Holstein-X-Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Holstein-Cross.deepv-snps.fst &

vcftools --gzvcf ../holPub.deepv-snps.filter.vcf.gz \
    --weir-fst-pop Holstein-X-Jersey \
    --weir-fst-pop Jersey \
    --fst-window-size 10000 \
    --fst-window-step 5000 \
    --maf 0.01 \
    --out Cross-Jersey.deepv-snps.fst &


cat ../holPub.samples.group | awk '{if($2 ~ "-"){$2="Cross"};print $2"\t"$1}' > holPub.id


for brd in Holstein Cross Jersey; do 
    cat holPub.id | awk -v Breed=$brd '($1 ~ Breed){print Breed,$2,Breed}' > fst.id.$brd
done

vars="pangenie-sv"

plink2 --threads $nthreads --vcf ../holPub.$var.filter.vcf.gz \
    --chr-set 29 --chr 1-29 --maf 0.01 --make-bed \
    --out holPub.$var.bfile

mv holPub.$var.bfile.fam holPub.$var.bfile.fam.bak
awk 'NR==FNR{a[$2]=$1}NR>FNR{$1=a[$2];print $0}' holPub.id holPub.$var.bfile.fam.bak > holPub.$var.bfile.fam

brd1=(Holstein Cross Jersey)
for i in {0..2};do
for j in {0..2};do
    if [[ $j > $i ]];
        then
        cat fst.id.${brd1[$i]} fst.id.${brd1[$j]} > fst.id.${brd1[$i]}.vs.${brd1[$j]}
        plink --bfile holPub.$var.bfile --fst --within fst.id.${brd1[$i]}.vs.${brd1[$j]} \
            --chr-set 29 --chr 1-29 --out fst.$var.${brd1[$i]}.vs.${brd1[$j]} &
        fi
    done
done

pangenie-var

vars="deepv-snps pangenie-snps"

for var in $vars; do
plink2 --threads $nthreads --vcf ../holPub.$var.filter.vcf.gz \
    --chr-set 29 --chr 1-29 --maf 0.01 --make-bed \
    --out holPub.$var.bfile

mv holPub.$var.bfile.fam holPub.$var.bfile.fam.bak
awk 'NR==FNR{a[$2]=$1}NR>FNR{$1=a[$2];print $0}' holPub.id holPub.$var.bfile.fam.bak > holPub.$var.bfile.fam

brd1=(Holstein Cross Jersey)
    for i in {0..2};do
    for j in {0..2};do
        if [[ $j > $i ]];
            then
           # cat fst.id.${brd1[$i]} fst.id.${brd1[$j]} > fst.id.${brd1[$i]}.vs.${brd1[$j]}
            plink --bfile holPub.$var.bfile --fst --within fst.id.${brd1[$i]}.vs.${brd1[$j]} \
                --window 10000 --window-step 5000 \
                --chr-set 29 --chr 1-29 --out fst.$var.${brd1[$i]}.vs.${brd1[$j]} &
            fi
        done
    done
done

