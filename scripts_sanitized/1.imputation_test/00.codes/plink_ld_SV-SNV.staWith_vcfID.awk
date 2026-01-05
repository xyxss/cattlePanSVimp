NR > 1{
  split($3,a,/:|-|_/);
  split($6,b,/:|-|_/);
    if((a[5] > 50 && b[5] <= 50) || (b[5] > 50 && a[5] <= 50)){
      if(a[5] > 50){
        if($5 >= a[2] && $5 <= a[3]){dis=0} else
        if(a[2] > $5 ){dis=a[2] - $5} else
        if($5 > a[3]){dis=$5 - a[2]} 
      } else {
        if($2 >= b[2] && $2 <= b[3]){dis=0} else
        if(b[2] > $2 ){dis=b[2] - $2} else
        if($2 > b[3]){dis=$2 - b[2]}  
      }
      SV1=a[4];
      SV2="SNP"
    } else
    if( (a[5] > 50) || (b[5] > 50) ){
      if(a[2] < b[2] && a[2]+a[5] < b[2]){dis=b[2] - (a[2]+a[5])} else
      if(a[2] > b[2]+b[5] && a[2]+a[5] > b[2]+b[5]){dis=a[2] - (b[2]+b[5])} else
      {dis=0}
      SV1=a[4];
      SV2=b[4];
      } else
    {
      if($5 > $2){
        dis=$5-$2
      } else {
        dis=$2-$5
      }
      SV1="SNP";
      SV2="SNP"
    }
    count[SV1"xx"SV2,int(dis/1000)]++
    sumR2[SV1"xx"SV2,int(dis/1000)]+=$7
  }
  END {
    printf "type\tdistance\tR2\n";

    for (subscript in count) {
        split(subscript, a, SUBSEP);
        #sumR2[SV1"xx"SV2,int(dis/1000)]

        printf "%s\t%s\t%s\n", a[1], a[2], sumR2[a[1], a[2]]/count[a[1], a[2]]
    }
  }