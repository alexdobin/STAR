# usage: awk -f calcUMIperCell.awk raw/matrix.mtx raw/barcodes.tsv filtered/barcodes.tsv | sort -k1,1rn > UMIperCell.txt
# output: column1 = total UMIs per cell 
#         column2 = 1 for cell that passed filtering, 0 otherwise

BEGIN {
    OFS="\t";
}

{
if (ARGIND==1) {
    if (FNR<4)
        next; #skip header

    umiCount[$2]+=$3;
     
} else if (ARGIND==2) {
    rawCB[$1]=FNR;
} else if (ARGIND==3) {
    filtCB[rawCB[$1]]=FNR;
}

}

END {
    for (ii in umiCount)
         print umiCount[ii], (ii in filtCB);
}
