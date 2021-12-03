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
    CB[$2][ARGIND]=1;
}

END {
    for (gc in C) {
        split(gc,a,"_");
        if (CB[a[2]][1]!=1 || CB[a[2]][2]!=1)
            continue; # only cells present in both 1 and 2 are considered

        n++;
        c2=C[gc][2]+0;
        c1=C[gc][1]+0;
        maxc12 = (c2>c1 ? c2:c1);
        if (maxc12==0)
            maxc12=1e-9;
        rd=(c2-c1)/maxc12;

        print c1, c2 > "counts.txt";
 
        MARD += (rd>0 ? rd:-rd);
        if (rd==1) NrdP1++;
        if (rd==-1) NrdM1++;
        if (rd==0) Nrd0++;

        cc[1]=log(c1+1); 
        cc[2]=log(c2+1); 
        sum12 += cc[1]*cc[2];
        for (ii=1; ii<=2; ii++) {
            sum1[ii] += cc[ii];
            sum2[ii] += cc[ii]^2;
        };
        
    };

    print "MARD", MARD/n, NrdM1, Nrd0, NrdP1, n-NrdM1-Nrd0-NrdP1;
    print "R2", (n*sum12-sum1[1]*sum1[2])^2 / (n*sum2[1]-(sum1[1])^2) / (n*sum2[2]-(sum1[2])^2);

};
