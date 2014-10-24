function cigarGenomicDist(cig)
{
        n=split(cig,L,/[A-Z]/)-1;
        split(cig,C,/[0-9]*/);
        g=0;
        for (ii=1;ii<=n;ii++) {//scan through CIGAR operations
                if (C[ii+1]!="S" && C[ii+1]!="I") { 
                        g+=L[ii];
                };
        };
        return g;
};
BEGIN {
        endTol=5;
};
{
if ( $7>=0 && $1==$4 && $3==$6 && (($3=="-" && $5>$2 && $5-$2<1000000) || ($3=="+" && $2>$5 && $2-$5<1000000)) )
{
	#print $1,$2,$5,$3,$7,$8,$9;
	#print $11,$11+cigarGenomicDist($12),$13,$13+cigarGenomicDist($14);
        if ( ($3=="+" && $11+endTol>$5 && $13+cigarGenomicDist($14)-endTol<=$2) \
          || ($3=="-" && $13+endTol>$2 && $11+cigarGenomicDist($12)-endTol<=$5) ) {  
               print $1,($3=="+"?$5:$2),($3=="+"?$2:$5),($3=="+"?"-":"+"),($7==0?0:3-$7),$8,$9;
	};
};
};
