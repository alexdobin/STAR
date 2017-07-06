BEGIN {
  while (getline < TrLengthFile) 
  {
      T[$1]=$2+0;
  } 
} 

{ 
  if (NR%2==1) 
  {
    s=$4;
    if (s==0)
      next;
    i=$9;
    l=T[$3];
    n=substr($12,6);

    #print s,i,l,n > "alignTr.txt"; 

    iN[i]+=1/n; 

    t1=int((s-1)/l*100)+1;t2=int((s+i-2)/l*100)+1; 
    for (ii=t1;ii<=t2;ii++) 
      C[ii]+=1/n;
   
    if (n==1)
    {
      for (ii=t1;ii<=t2;ii++) 
        Cu[ii]+=1;

      iNu[i]++;
    }
  } 
}
END {
  for (ii in iN) 
    print ii,iN[ii] > "insertHist.txt"; 
  for (ii in iNu) 
    print ii,iNu[ii] > "insertHistUnique.txt"; 
  for (ii=1;ii<=101;ii++) 
    print ii,C[ii]+0,Cu[ii]+0 > "coverage.txt"
}
