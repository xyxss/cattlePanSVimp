BEGIN{
FS=OFS="\t"
}
{
    if($1 ~ "#"){
        print
    } else {
        for(i=10;i<=NF;i++){
            if($i=="."){$i="./."}
        }
        xi++;
        split($8, a, ";");
        for(i in a){
            if(a[i] ~ /^SVTYPE/){
                split(a[i], b, "=");
                typ=b[2]
            }
            if(a[i] ~ /^SVLEN/){
                split(a[i], b, "=");
                leng=b[2]
            }
            if(a[i] ~ /^ID/){
                split(a[i], b, "=");
                id=b[2]
            }
        }
        if(leng < 50){
            typ="s"typ
        }
        $3=$1":"$2":"PG"-"typ"_"leng"_"xi":"id
        print 
    }
}
