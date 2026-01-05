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
    typ="SNP";
    leng=1;
    pos=$2;
    len1=length($4);
    len2=length($5);
    if(len2 > 1){
      str1=$4;
      str2=$5;
      len = (len1 < len2) ? len1 : len2;
      for (i = 1; i <= len; i++) {
        if (substr(str1, i, 1) != substr(str2, i, 1)) {
          $4=substr(str1, i, 1) 
          $5=substr(str2, i, 1);
          $2=$2+i-1;
          break;
        }
      }
    }
    $3=$1":"$2":"$2"-"typ"_"leng"_"NR":"pos
    print
}