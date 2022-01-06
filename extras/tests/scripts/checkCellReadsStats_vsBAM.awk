# usage awk -f checkCellReadsStats.awk <(samtools view Aligned.sortedByCoordinate.bam)

BEGIN {
    OFS="\t";
}

(ARGIND==1) {
    if ($3=="*") {#unmapped
        CB=substr($(NF-1),6);
        if (CB!="-")
            R[$1]["CB"]=CB;
        next;
    };


    CB=substr($26,6);
    if (CB=="-")
        next;
    GX=substr($21,6);
    gx=substr($23,6);

    #print $1,CB,GX,gx;
    R[$1]["mapped"]++;

    R[$1]["CB"]=CB;
    if (GX!="-")
        R[$1]["GX"]++;    
    if (gx!="-")
        R[$1]["gx"]++;    
    if ($3=="chrM")
        R[$1]["chrM"]++;
}

(ARGIND==2) {
    if (FNR==1) {
        for (ff=1; ff<=NF; ff++) {
            fn[ff]=$ff;
            fi[$ff]=ff;
        };
        next;
    };

    for (ff=2; ff<=NF; ff++) {
        crs[$1][fn[ff]]=$ff;
    };
}

END {
    for (r in R) {
        CB=R[r]["CB"];
        c[CB]["cbMatch"]++;
        if (R[r]["GX"]>0) {
            c[CB]["featureU"]++;
        } else if (R[r]["gx"]>0) {
            c[CB]["featureM"]++;
        };

        if (R[r]["mapped"]==1)
            c[CB]["genomeU"]++;
        if (R[r]["mapped"]>1)
            c[CB]["genomeM"]++;
        if (R[r]["chrM"]>0)
            c[CB]["mito"]++;

    };

    split("genomeU genomeM featureU featureM mito", fnames);
    for (CB in c) {
        #print CB,c[CB]["cbMatch"],crs[CB]["cbMatch"];
        for (f=1; f<=length(fnames); f++) {
             ff=fnames[f];
             if (crs[CB][ff]!=c[CB][ff])
                 print CB, ff, crs[CB][ff]+0, c[CB][ff]+0;
        };
    };


    # write out the file - not needed
    exit;
    split("cbMatch genomeU genomeM featureU featureM mito", fnames);
    for (CB in c) {
        printf CB > "CellReads.txt";
        for (f=1; f<=length(fnames); f++) {
             ff=fnames[f];
             printf "\t" c[CB][ff]+0 > "CellReads.txt";
        };
        printf "\n" > "CellReads.txt";
    };
}
