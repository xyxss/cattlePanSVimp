BEGIN {
    FS=OFS="\t"
}
$1 ~ /#/{
    print
    next
}
{
    for(i=10;i<=NF;i++){
        if($i=="."){$i="./."}
    }
    print
}