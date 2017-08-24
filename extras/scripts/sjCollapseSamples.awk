# usage:
# awk -f sjCollapseSamples.awk /path/to/all/*/SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > SJ.all
# output columns:
# 1-6 - same as in SJ.out.tab
# 7   - total number of unique mappers
# 8   - total number of multi-mappers
# 9   - max overhang
# 10  - number of samples the junction was detected in

BEGIN {
    OFS="\t";
}

{
    sj=$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6; 
    nSamples[sj]++;
    nU[sj]+=$7;
    nM[sj]+=$8;
    if (nO[sj]<$9) nO[sj]=$9;
};

END {
    for (sj in nSamples) print sj,nU[sj],nM[sj],nO[sj],nSamples[sj];
}

