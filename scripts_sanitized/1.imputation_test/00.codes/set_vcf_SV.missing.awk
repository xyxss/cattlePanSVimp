BEGIN {
    FS=OFS="\t"
}
$1 ~ /##/{
    print
    next
}
$1 ~ /#/{
    for(i=10;i<=NF;i++){
        sample_ids[i] = $i
    }
    print
    next
}   
{   
    split($3,a,/:|-|_/);
        if(a[5] > 50){
            for(i=10;i<=NF;i++){
                print $3,$i > "set.missing"rate"/set."sample_ids[i]".missing"rate".IDs"
                $i="./."
            }
        };
    print
}