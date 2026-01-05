



mkdir -p /90daydata/bull_age/liu.yang/imputation/imputation
cd /90daydata/bull_age/liu.yang/imputation/imputation

work_path="/90daydata/bull_age/liu.yang/imputation/imputation"

#### sofrware
cd /home/liu.yang/software
conda create -n imp -c conda-forge -c bioconda plink plink2 eagle bcftools vcftools -y

wget http://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar

java -jar beagle.06Aug24.a91.jar gt=test.06Aug24.a91.vcf.gz out=out.gt
java -jar beagle.06Aug24.a91.jar ref=ref.06Aug24.a91.vcf.gz gt=target.06Aug24.a91.vcf.gz out=out.ref

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
#wget https://github.com/statgen/Minimac4/releases/download/v4.1.2/minimac4-4.1.2-Linux.sh

wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh


ref_path="/90daydata/bull_age/liu.yang/ref/"
fa="ARS_UCD_v2.0.fa"
ref_fa=$ref_path/$fa

minimac4=/home/liu.yang/software/minimac4-4.1.6-Linux-x86_64/bin/minimac4
beagle=/home/liu.yang/software/beagle.06Aug24.a91.jar
eagle=/home/liu.yang/software/Eagle_v2.4.1/eagle

$minimac4 --compress-reference reference.{sav,bcf,vcf.gz} > reference.msav
$minimac4 reference.msav target.vcf.gz -o imputed.vcf.gz


input='/home/yangliu/3.imp/pig_SV_1060.imputation/pig_SV_1060.chr1.vcf.gz'
geneticmap='/home/yangliu/3.imp/pig_SV_1060.imputation/Sus_scrofa.Sscrofa11.1.dna.toplevel.nochr.fa.map.txt'
chromosome="18"
multiple_SV="Ture"


input_pre_bfile="/home/yangliu/3.imp/SV.all.bfile"
input_pre_vcf_gz="/home/yangliu/3.imp/pig_sv.in.vcf.gz"
plink_keep="/home/yangliu/3.imp/plink.keep"
pre="pig_SV_1060"
chr="10"
chrs=18


plink --bfile $input_pre_bfile --keep $plink_keep --const-fid 0 --chr 1-$chrs --recode vcf-iid --out $pre.chr$chrs


generate_geneticmap_1cM_Mb() {
  cat $ref_fa.fai |
    awk -v chrs=$1 '$1 <= chrs{
    max=int($2/1e6);
    count=0;
    while(count < max) {
      print $1,int($2/max*count),1,count
        count++;
  }}' | sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \
    >${ref_fa//*\//}.map.txt
}

mkdir -p /home/yangliu/3.imp/$pre.imputation
cd /home/yangliu/3.imp/$pre.imputation

if [[ ! -f ${ref_fa//*\//}.map.txt ]]; then
  generate_geneticmap_1cM_Mb $chrs
fi

chrss="1
10
18"

tabix -p vcf $input_pre_vcf_gz

bcftools view $input_pre_vcf_gz --regions $chr |
    awk 'BEGIN{FS=OFS="\t"}{if($1 !~ "#"){
    a[$1,$2]++;
    if(a[$1,$2] > 1){
      i=1
      while(a[$1,$2 + i] > 1){
        i++
      }
      $2+=i
      a[$1,$2]++
    }
  }
  print}' | bgzip -c >$pre.chr$chr.vcf.gz

  cat ${ref_fa//*\//}.map.txt | awk -v chr=$chr '{if(NR == 1 || $1 == chr){print}}' \
    >${ref_fa//*\//}.$chr.map.txt


  bcftools view -i 'F_MISSING < 0.1 && MAF > 0.01' $pre.chr$chr.vcf.gz -O z -o $pre.chr$chr.vcf.filter.maf0.01_fmiss0.01.vcf.gz

  bcftools view -i 'F_MISSING < 0.1 && MAF > 0.05' $pre.chr$chr.vcf.gz -O z -o $pre.chr$chr.vcf.filter.maf0.05_fmiss0.05.vcf.gz

  bcftools view -i 'F_MISSING < 0.1 && MAF > 0.2' $pre.chr$chr.vcf.gz -O z -o $pre.chr$chr.vcf.filter.maf0.05_fmiss0.2.vcf.gz
