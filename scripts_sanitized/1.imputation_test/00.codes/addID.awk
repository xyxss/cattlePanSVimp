BEGIN{
FS=OFS="\t"
}
$1 ~ /##/{
    print
    next
}
$1 ~ /#/{
    print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    print "##INFO=<ID=.,Number=1,Type=String,Description=\"aaa\">"
    print
    next
}
{
    xi++;
    altlen=length($5);
    reflen=length($4);
    leng=0;
    if(altlen > reflen && reflen == 1){
        typ="INS"
        leng=altlen-1
    } else if(altlen < reflen && altlen == 1){
        typ="DEL"
        leng=reflen-1
    } else if(altlen == reflen && reflen == 1){
        typ="SNP"
        leng=1
    } else {
        typ="COMPLEX"
        leng=altlen-1
    } 
    
    if(leng < 50){
        typ="s"typ
    }
    
    $3=$1":"$2":"PG"-"typ"_"leng"_"xi":"$3
    print
}