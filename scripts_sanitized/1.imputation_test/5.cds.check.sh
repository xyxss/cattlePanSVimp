

cd /90daydata/bull_age/liu.yang/imputation/imputation/6.cds.check



vcf=/90daydata/bull_age/liu.yang/imputation/imputation/1.holPub_imp/run_files/Holstein750_Thin0k_001.pangenie-snps.chr10.vcf.gz
ref_rm=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_repeat
ref_gff=/90daydata/bull_age/liu.yang/ref/ARS_UCD_v2.0.ref_gff


rms="None RM DNA LTR Low_complexity LINE srpRNA rRNA Unknown RNA RC scRNA SINE Satellite Simple_repeat tRNA snRNA"
gfs="nongene
intron
antisense_RNA
CDS
C_gene_segment
D_loop
exon
gene
lnc_RNA
miRNA
mRNA
ncRNA
primary_transcript
pseudogene
RNase_MRP_RNA
rRNA
snoRNA
snRNA
SRP_RNA
telomerase_RNA
transcript
tRNA
V_gene_segment"

cat > hetRate.awk << 'EOF'
{
    total_sites++;  
    for (i=3; i<=NF; i++) {
        if ($i ~ /[0-9]\/[0-9]/ && substr($i, 1, 1) != substr($i, 3, 1)) {
            heterozygous++;
        }
    }
} END {
    heterozygosity_rate = total_sites > 0 ? heterozygous / total_sites : 0;
    print "Total sites:", total_sites;
    print "Heterozygous sites:", heterozygous;
    printf "Heterozygosity rate: %.4f\n", heterozygosity_rate;
}
EOF

bcftools stats $vcf > stats.all.txt

rm stats.rm.txt

for rm in $rms; do
bcftools view -R $ref_rm/rm.$rm.bed $vcf | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | \
    awk -f hetRate.awk | awk -v rm=$rm '{a[NR]=$3};END{print rm"\t"a[1],a[2],a[3]}' >> stats.rm.txt &
done

rm stats.gf.txt 
bcftools view $vcf | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | \
    awk -f hetRate.awk | awk -v rm="All" '{a[NR]=$3};END{print rm"\t"a[1],a[2],a[3]}' >> stats.gf.txt &
for gf in $gfs; do
bcftools view -R $ref_gff/gff.$gf.bed $vcf | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | \
    awk -f hetRate.awk | awk -v rm=$gf '{a[NR]=$3};END{print rm"\t"a[1],a[2],a[3]}' >> stats.gf.txt &
done

cat stats.all.txt | awk '$1 == "AF" {print "All",$3,$4}' > afc.gf.txt 
for gf in $gfs; do
bcftools stats -R $ref_gff/gff.$gf.bed $vcf | awk '$1 == "AF" {print rm,$3,$4}' rm=$gf >> afc.gf.txt &
done



