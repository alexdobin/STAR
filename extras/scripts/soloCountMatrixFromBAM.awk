# usage: 
# samtools view Aligned.sortedByCoord.out.bam | awk -v fileWL=Solo.out/Gene/raw/barcodes.tsv -v fileGenes=Solo.out/Gene/raw/features.tsv \
# soloCountMatrixFromBAM.awk | sort -k2,2n -k1,1n > matrix.mtx


function getTag(tag)
{
    tagOut=$0;
    if (gsub(".*" tag,"",tagOut)==0)
        return 0;
    gsub("\t.*","",tagOut);
    return tagOut;
}

BEGIN {

    stderr="/dev/stderr";

    ii=0;
    while (getline < fileWL) {
        ii++;
        WL[$1]=ii;
    };
    #print length(WL) > stderr;

    #getline < fileGenes; #skip header
    ii=0;
    while (getline < fileGenes) {
        ii++;
        geneID[$1]=ii;
    };
    #print length(geneID) > stderr;
}

{
    #if ((NR%1000000)==0) printf NR/1000000 " "$3 " " nTot+0 "  "  > "/dev/stderr";

    GX=getTag("GX:Z:");
    if (GX=="0" || GX=="-")
        next;

    CB=getTag("CB:Z:");
    if (CB=="0" || CB=="-")
        next;

    UB=getTag("UB:Z:");
    if (UB=="0" || UB=="-")
        next;

    # this is needed for CR
    #if (substr(CB,length(CB)-1,1)=="-")
    #    CB=substr(CB,1, length(CB)-2);


    cb=WL[CB];
    ge=geneID[GX];


    nTot++;
    U[cb][ge][UB]++;
      
}

END {
    for (cb in U) {
        for (ge in U[cb]) {
            print ge,cb,length(U[cb][ge]);
        }
    }   
}
