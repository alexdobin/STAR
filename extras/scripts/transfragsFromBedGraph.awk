BEGIN {
    OFS="\t";
    getline;
    start1=$2;
    end1=$3;
    chr1=$1;
    s=$4;
} 

{
    if ($2!=end1 && start1>0) {
        printf "%s\t%i\t%i\t%20.5f\n", chr1,start1,end1,s; 
        chr1=$1; 
        start1=$2; 
        end1=$3; 
        s=$4
    } else {
        s+=$4;
        end1=$3; 
        #print end1,start1,$2
    };
}

END {
    printf "%s\t%i\t%i\t%20.5f\n", chr1,start1,end1,s;
}
