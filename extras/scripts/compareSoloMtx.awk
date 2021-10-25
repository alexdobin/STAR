BEGIN {
    OFS="\t";
}

BEGINFILE {
    skipLines=0;
}

{
    if ($1~/^%/) next;
    if (skipLines==0) {
        skipLines=1;
        next;
    };
    C[$1 "_" $2][ARGIND]=$3;
}

END {
    for (gc in C) {
        c2=C[gc][2]+0;
        c1=C[gc][1]+0;
        rd=(c2-c1)/(c2>c1 ? c2:c1);
 
        MARD += (rd>0 ? rd:-rd);
        if (rd==1) NrdP1++;
        if (rd==-1) NrdM1++;
        if (rd==0) Nrd0++;
     
        sum12 += c1*c2;
        for (ii=1; ii<=2; ii++) {
            sum1[ii] += C[gc][ii];
            sum2[ii] += C[gc][ii]^2;
        };
        
    };

    n=length(C);
    print "MARD", MARD/n, NrdM1, Nrd0, NrdP1;
    print "R2", (n*sum12-sum1[1]*sum1[2])^2 / (n*sum2[1]-(sum1[1])^2) / (n*sum2[2]-(sum1[2])^2);

};
