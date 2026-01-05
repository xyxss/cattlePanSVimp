cd /90daydata/bull_age/liu.yang/imputation/imputation/4.hd_imp

wget https://webdata.illumina.com/downloads/productfiles/bovinehd/bovinehd-b1-annotation-file.zip
conda install bioconda::ucsc-liftover
wget https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates/UMC_marker_names_180910.zip

wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver

chmod +x ./filePath/utility_name
./filePath/utility_name

https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/liftOver/

#ref https://groups.google.com/a/soe.ucsc.edu/g/genome/c/pmpGHtsyiP4


Here are the details for the two assemblies:
GCF_002263795.3_ARS-UCD2.0 cattle (Hereford L1 Dominette 01449 42190680 v2.0 2023 USDA)
GCF_000003055.6_Bos_taurus_UMD_3.1.1 cattle (Hereford UMD_3.1.1 2014)

The LiftOver files are available here:
https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/liftOver/GCF_002263795.3ToGCF_000003055.6.over.chain.gz
https://hgdownload.soe.ucsc.edu/hubs/GCF/000/003/055/GCF_000003055.6/liftOver/GCF_000003055.6ToGCF_002263795.3.over.chain.gz

And you can view these genomes and the chain tracks (that display these alignments) here:
https://genome.ucsc.edu/h/GCF_002263795.3
https://genome.ucsc.edu/h/GCF_000003055.6

I hope this is helpful. If you have any further questions, please reply to gen...@soe.ucsc.edu. All messages sent to that address are archived on a publicly-accessible Google Groups forum. If your question includes sensitive data, you may send it instead to genom...@soe.ucsc.edu.

wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/003/055/GCF_000003055.6/liftOver/GCF_000003055.6ToGCF_002263795.3.over.chain.gz

#cat BovineHD_B1.csv | awk -F, '$9 == "UMD_3.1" && $10 > 0 {print "chr"$10"\t"$11"\t"$11}' > BovineHD.UMD_3.1.bed

cd /90daydata/bull_age/liu.yang/imputation/hd_gwas

# liftover =
{
cat BovineHD_B1.csv | awk -F, '$9 == "UMD_3.1" {print $10"\t"$11"\t"$2}' > BovineHD.UMD_3.1.pos

awk 'FNR==NR{a[$1]=$2;next} {print a[$1]"\t"$2"\t"$2"\t"$3}' /90daydata/bull_age/liu.yang/ref/GCF_000003055.6_UMD3.1.chr.tr BovineHD.UMD_3.1.pos  | awk '$1 ~ "AC"' > BovineHD.GCF_000003055.6_UMD3.1.bed
#775032

liftOver BovineHD.GCF_000003055.6_UMD3.1.bed GCF_000003055.6ToGCF_002263795.3.over.chain GCF_002263795.3_ARS-UCD2.0.bed unMapped.bed

#ARS-UCD2.0.bed
#774316

awk 'FNR==NR{a[$2]=$1;next} {print a[$1]"\t"$2"\t"$2"\t"$4}' /90daydata/bull_age/liu.yang/ref/GCF_002263795.3_ARS-UCD2.0.chr.tr GCF_002263795.3_ARS-UCD2.0.bed  > ARS-UCD2.0.bed

cat GCF_002263795.3_ARS-UCD2.0.bed | awk '{print $1"\t"$2}' > GCF_002263795.3_ARS-UCD2.0.pos
}

bcftools view --threads 4 -v snps holPub.pangenie.vcf.gz \
    -Oz -o holPub.pangenie-snps.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-snps.vcf.gz

## hdsnps
bcftools view --threads $nthreads -T GCF_002263795.3_ARS-UCD2.0.pos \
    ../1.holPub_imp/holPub.pangenie-snps.vcf.gz \
    -Oz -o holPub.pangenie-hdsnps.vcf.gz && 
tabix -f -p vcf holPub.pangenie-hdsnps.vcf.gz

work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp
run_file=$work_dir/run_files
code=$work_dir/00.codes

## svhdsnps
bcftools concat --threads $nthreads \
    --allow-overlaps \
    ../1.holPub_imp/holPub.pangenie-sv.vcf.gz \
    holPub.pangenie-hdsnps.vcf.gz |
    awk -f $code/haploid2diploid.awk |
    bgzip -c > holPub.pangenie-svhdsnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-svhdsnps.vcf.gz


work_dir=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp
run_file=$work_dir/run_files
code=$work_dir/00.codes
mkdir -p $run_file

mkdir -p $work_dir/imp_runs
chr=10

bcftools view --threads $nthreads holPub.pangenie-hdsnps.vcf.gz \
    -r $chr -S ../1.holPub_imp/Holstein750 --force-samples | 
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-hdsnps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-hdsnps.chr$chr.vcf.gz

bcftools view --threads $nthreads holPub.pangenie-svhdsnps.vcf.gz \
    -r $chr -S ../1.holPub_imp/Holstein750 --force-samples |
    awk -f $code/haploid2diploid.awk |
    bgzip -c -@ $nthreads > $run_file/Holstein750_Thin0k.pangenie-svhdsnps.chr$chr.vcf.gz &&
    tabix -f -p vcf -@ $nthreads $run_file/Holstein750_Thin0k.pangenie-svhdsnps.chr$chr.vcf.gz

chr=10
nt=2
proc=imp_all.sh
for sub in pangenie-svhdsnps pangenie-hdsnps; do
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

## xxxxx

function grp_vcf () {
    local vcf_file=$1
    local out_file=$2
    for grp in Holstein250 HolsteinRelated250 Jersey250 MultiBreed750; do
        bcftools view --threads $nthreads -r $chr \
             -S ../1.holPub_imp/$grp $vcf_file |
             awk -f $code/haploid2diploid.awk |
            bgzip -c ${out_file/holPub/$grp} && 
        tabix -f -p vcf ${out_file/holPub/$grp}
        echo ${out_file/holPub/$grp}
    done
}

grp_vcf holPub.pangenie-hdsnps.vcf.gz run_files/holPub_Thin0k.pangenie-hdsnps.chr$chr.vcf.gz
grp_vcf holPub.pangenie-svhdsnps.vcf.gz run_files/holPub_Thin0k.pangenie-svhdsnps.chr$chr.vcf.gz


chr=10
nt=2
proc=imp_all.sh
for sub in pangenie-hdsnps pangenie-svhdsnps; do
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

