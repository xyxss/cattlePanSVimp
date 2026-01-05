
/90daydata/bull_age/liu.yang/imputation/pangenie_WGS/cdcb_bam/


bcftools view --threads 4 -v indels \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.indel.awk |
    bgzip -c > holPub.pangenie-indels.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-indels.vcf.gz &

bcftools view --threads 4 -v other,mnps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > holPub.pangenie-complex.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-complex.vcf.gz &

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-complex.vcf.gz \
    holPub.pangenie-indels.vcf.gz \
    -Oz -o holPub.pangenie-nonsnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-nonsnps.vcf.gz

zcat holPub.pangenie-nonsnps.vcf.gz | 
    awk '$1 ~ /#/{print; next}{split($3,a,/:|-|_/);if(a[5] > 50) {print}}' |
    bgzip -c > holPub.pangenie-sv.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-sv.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-sv.vcf.gz -Oz -o holPub.pangenie-sv.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-sv.filter.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-snps.vcf.gz -Oz -o holPub.pangenie-snps.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-snps.filter.vcf.gz

bcftools view --threads 4 -v snps \
    -r $(seq 1 29 | tr '\n' ',') -i 'AC > 0' \
    holPub.pangenie.vcf.gz |
    awk -f  $code/pan.snp.awk |
    bgzip -c > holPub.pangenie-snps.vcf.gz && 
    tabix -f -p vcf holPub.pangenie-snps.vcf.gz

bcftools concat --threads $nthreads \
    --allow-overlaps \
    holPub.pangenie-nonsnps.vcf.gz \
    holPub.pangenie-snps.vcf.gz \
    -Oz -o holPub.pangenie-var.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-var.vcf.gz

bcftools view --threads 4 -i 'F_MISSING < 0.1 && MAF > 0.01' \
    holPub.pangenie-var.vcf.gz -Oz -o holPub.pangenie-var.filter.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-var.filter.vcf.gz

bcftools view --threads 4 -v mnps -r $(seq 1 29 | tr '\n' ',') holPub.pangenie.vcf.gz | 
    awk -f  $code/pan.complex.awk |
    bgzip -c > holPub.pangenie-mnps.vcf.gz &&
    tabix -f -p vcf holPub.pangenie-mnps.vcf.gz &

