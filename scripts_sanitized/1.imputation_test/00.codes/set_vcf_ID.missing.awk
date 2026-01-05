BEGIN {
    FS=OFS="\t"
    srand(9)
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
    for(i=10;i<=NF;i++){
        if($i=="."){$i="./."}
        if(rand()<rate && $i!="./.") {
            print $3,$i > "set.missing"rate"/set."sample_ids[i]".missing"rate".IDs"
            $i="./."
        }
    }
    print
}