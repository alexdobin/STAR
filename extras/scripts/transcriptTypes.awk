# requires "trTypes.txt" - a file with transcript types
# e.g. for Gencode GTF
# awk '$3=="transcript" {a=$0; gsub(/.*transcript_id "/,"",a);gsub(/".*/,"",a);;b=$0; gsub(/.*gene_type "/,"",b);gsub(/".*/,"",b); print a,b}' Gencode.gtf > trTypes.txt

BEGIN {
  while (getline < "trTypes.txt") {
    tT[$1]=$2;
  };
  OFS="\t";
  rt[1]=0; #declare array
  delete rt;
} 

{
  if ($1!=r) {
    if (length(rt)==1) {#only if read overlaps one trType
      for (tt in rt) {
        if (tt=="") print r;
        nT[tt]++;
      };
    };
    delete rt;
    r=$1;
  };

  if ($3 in tT) {    
    rt[tT[$3]]=1;
  };
};
 
END {
  for (tt in nT) {
    print tt, nT[tt];
  };
};

