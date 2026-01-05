
cd /90daydata/bull_age/liu.yang/imputation/sbcp_imp/genotypes


vcf=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/holPub.pangenie-var.vcf.gz

bcftools view --threads $nthreads $vcf \
    -S Holsteinandrelated.samples --force-samples \
    -Oz -o Holsteinandrelated.pangenie-var.vcf.gz &&
    tabix -f -p vcf -@ $nthreads Holsteinandrelated.pangenie-var.vcf.gz


zcat Holsteinandrelated.pangenie-var.vcf.gz | awk '$1 !~ "#" {split($3,a,/[:\-_]/);print a[1]"\t"a[2]"\t"a[3]"\t"$3}' > Holsteinandrelated.pangenie-var.pos
15727951
cat Holsteinandrelated.pangenie-var.pos | awk '{split($4,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}'
SNP 12512405
DEL 19770
sDEL 1250141
INS 37684
sINS 1318751
COMPLEX 35605
sCOMPLEX 553595

bedtools intersect -a  Holsteinandrelated.pangenie-var.pos \
    -b /90daydata/bull_age/liu.yang/ref/GCF_002263795.3_ARS-UCD2.0_genomic_gaps.autochr.txt \
    -r -v -f 0.1 -wa \
    >  Holsteinandrelated.pangenie-var.nogap.pos.bed
15727951

bedtools intersect -a  Holsteinandrelated.pangenie-var.pos \
    -b /90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat/rm.Satellite.bed \
    -r -v -f 0.1 -wa \
    >  Holsteinandrelated.pangenie-var.noSate.pos.bed
cat Holsteinandrelated.pangenie-var.noSate.pos.bed | awk '{split($4,a,/[:\-_]/);b[a[4]]++} END {for(i in b) {print i,b[i]}}'
15727718


seq 1 29 | while read chr; do
    sbatch -A bull_age \
        --job-name=chr$chr.refIm \
        --cpus-per-task=48 \
        --error=logs/refIm.chr$chr.err \
        --output=logs/refIm.chr$chr.out \
        --wrap="
        bash refImp.sh $chr
        "
done