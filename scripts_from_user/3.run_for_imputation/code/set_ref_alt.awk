BEGIN {
    FS=OFS="\t"
}
$1 ~ /#/{
    print
    next
}
{
    $4="A"
    $5="C"
    print
}