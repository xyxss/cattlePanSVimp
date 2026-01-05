NR > 1{
    if(a[5] > 50){
      id=$3
      idx=$6
      split($3,a,/:|-|_/);
      if($5 >= a[2] && $5 <= (a[2]+a[5])){dis=0} else
      if(a[2] > $5 ){dis=a[2] - $5} else
      if($5 > (a[2]+a[5])){dis=$5 - a[2]} 
    } else {
      id=$6
      idx=$3
      split($6,a,/:|-|_/);
      if($2 >= a[2] && $2 <= (a[2]+a[5])){dis=0} else
      if(a[2] > $2 ){dis=a[2] - $2} else
      if($2 > (a[2]+a[5])){dis=$2 - a[2]}  
    }
    SV1=a[4];
    x=int(($7*10+0.5))/10
    sdis=int(dis/1000+0.5)
    svt[x,SV1,id,sdis]++
    svt[x,"SNP",idx,sdis]++
    #print x,SV1,id,sdis
  } END {
    #printf "type\tdistance\tR2\n";
    for (subscript in svt) {
      split(subscript, a, SUBSEP);
      count[a[1],a[2],a[4]]++
    }

    for (subscript in count) {
      split(subscript, a, SUBSEP);
      printf "%s\t%s\t%s\t%s\n", a[2], (a[2]+a[5]), a[1], count[a[1],a[2],(a[2]+a[5])] 
    }
  }