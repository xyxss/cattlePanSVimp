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
    print;
    next
}
{

    altlen=length($5);
    reflen=length($4);
    leng=0;
    typ="COMPLEX"
    if(altlen > reflen){
        leng=altlen-1
    } else  {
        leng=reflen-1
    }
    #print $1,$2,$2+leng,leng,typ,"-",typ"-"NR 
    if(leng <= 50) {typ = "s"typ}
    $3=$1":"$2":"$2+leng"-"typ"_"leng"_"NR":"$3
    print
}