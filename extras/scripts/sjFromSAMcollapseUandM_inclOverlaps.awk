BEGIN {
   OFS="\t";
   mapqU=255;
}
{
if (substr($1,1,1)!="@") {

    m=and($2,0x80)/0x80+1;

    if ($1!=readNameOld) delete readSJs;
    readNameOld=$1;

    n=split($6,L,/[A-Z]/)-1;
    split($6,C,/[0-9]*/);
    t=1;g=$4;
    for (k=1;k<=n;k++) {#scan through CIGAR operations
        if (C[k+1]=="S" || C[k+1]=="I") {
           t+=L[k];
        } else if (C[k+1]=="D") {
           g+=L[k];
        } else if (C[k+1]=="N") {
           sj1=$3 "\t" g "\t" g+L[k]-1;
           readSJs[sj1]++;

           if (readSJs[sj1]==1) {#only count this junction if it has nto been counted for the same read
               SJ[sj1]=1;
               if ($5>=mapqU) {
                   SJu[sj1]++;
               } else {
                   SJm[sj1]++;
               };
           };

           if ($5>=mapqU) {
               SJu1[sj1]++;
           } else {
               SJm1[sj1]++;
           };

           g+=L[k];

        } else { # M operation
           g+=L[k];
           t+=L[k];
        };
    };
};
};
END {

for (ii in SJ) {
    print ii, SJu[ii]+0, SJm[ii]+0, SJu1[ii]+0, SJm1[ii]+0;
};

};
