# usage awk -f checkCellReadsStats.awk CellReads.stats barcodes.tsv matrix.mtx
function funAbs(x) {
    return (x>0 ? x : -x);
};

BEGIN {
    OFS="\t";
}

(ARGIND==1) {
if (NR==1) {
    for (ff=1; ff<=NF; ff++) {
        fn[ff]=$ff;
        fi[$ff]=ff;
    };
    next;
};

nuu=$(fi["nUMIunique"]);
ngu=$(fi["nGenesUnique"]);
num=$(fi["nUMImulti"]);
ngm=$(fi["nGenesMulti"]);

if (nuu>0) {
    nGenes[$1]=ngu;
    nUMI[$1]=nuu;
};
if (nuu+num>0) {
    nGenesTot[$1]=ngu+ngm;
    nUMItot[$1]=nuu+num;
};

}

(ARGIND==2) {#barcodes.tsv
    CB[FNR]=$1;
}

(ARGIND==3 && FNR>3) {#matrix mtx
    matGenes[CB[$2]]++;
    matUMI[CB[$2]]+=$3;
}

(ARGIND==4 && FNR>3) {#matrix mtx
    matGenesTot[CB[$2]]++;
    matUMItot[CB[$2]]+=$3;
}


END {

    #for (cb in matUMI)
    #    print cb,matUMI[cb],matGenes[cb];

    #exit;

    print length(nGenes), length(matGenes);

    for (cb in matGenes) {
        if (nGenes[cb]!=matGenes[cb])
            print "G",cb,nGenes[cb]+0,matGenes[cb]+0;
        if (nUMI[cb]!=matUMI[cb])
            print "U",cb,nUMI[cb]+0,matUMI[cb]+0;
    };

    if (ARGC>4) {
        print length(nGenesTot), length(matGenesTot);

        for (cb in matGenesTot) {
            if (nGenesTot[cb]!=matGenesTot[cb])
                print "G",cb,nGenesTot[cb]+0,matGenesTot[cb]+0;
            if (funAbs(nUMItot[cb]-matUMItot[cb])>0.0001)
                print "U",cb,nUMItot[cb]+0,matUMItot[cb]+0;
        };
    };

}
