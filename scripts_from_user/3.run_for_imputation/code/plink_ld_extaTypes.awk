NR==1{
    print $0" SVTYPE VARTYPE DISTANCE"; 
    next
}
{
    split($3,a,/:|-|_/)
    split($6,b,/:|-|_/)

    if(a[5] > 50 && b[5] > 50 ){
        typ=a[4]"-"b[4]
        svtyp="sv-sv"
        if(a[2] < b[2] && a[2]+a[5] < b[2]){
            dis=b[2] - (a[2]+a[5])
        } else 
        if(a[2] > b[2]+b[5] && a[2]+a[5] > b[2]+b[5]){
            dis=a[2] - (b[2]+b[5])
        } else {dis=0}
    } else if(a[5] > 50){
        typ=a[4]"-snv"
        svtyp="sv-snv"
        if(a[2] > b[2] ){dis=a[2] - b[2]} else 
        if(b[2] > (a[2]+a[5])){dis=b[2] - (a[2]+a[5])} else {dis=0}
    } else if(b[5] > 50 ){
        typ=b[4]"-snv"
        svtyp="sv-snv"
        if(b[2] > a[2] ){dis=b[2] - a[2]} else
        if(a[2] > (b[2]+b[5])){dis=a[2] - (b[2]+b[5])} else {dis=0}
    } else {
        typ="snv-snv"
        svtyp="snv-snv"
        if(a[2] > b[2]){dis=a[2] - b[2]} else {dis=b[2] - a[2]}
        }
    print $0,typ,svtyp,dis
}