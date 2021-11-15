# usage awk soloExtractFiltCells.awk barcodes.tsv barcodes_filtered.tsv matrix.mtx > matrix_filtered.mtx
BEGIN {
    nSkip=1; # number of cells to skip after comments
}

{
    if (ARGIND==1) {# read barcodes
        CBind[$1]=NR;
    } else if (ARGIND==2) {
        n++;
        CB[CBind[$1]]=n; # new indexes
    } else if ($1~/^%/){
        print;
        next;
    } else {
        if (nSkip>0) {
            nFeat=$1;
            nSkip--;
            next;
        }; 
        if ($2 in CB) {
            nLines++;
            outAll = outAll sprintf($1 " " CB[$2] " " $3 "\n");
        };
    };
}

END {
      print nFeat, length(CB), nLines;
      print outAll;
};
