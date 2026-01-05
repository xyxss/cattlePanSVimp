BEGIN {
    srand(1)
    } {
    max=NF-9;
    total=int((NF-9)/10);
    count=0;
    while(count < total) {
      rnd = int(rand() * max);
      if ( array[rnd] == 0 ) {
        count++;
        array[rnd]++;
      }
    }
    for ( i=1; i<=max; i++) {
      if ( array[i] )
        print $i;
    }
  }