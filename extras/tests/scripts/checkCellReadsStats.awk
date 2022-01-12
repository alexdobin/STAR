# usage awk -f checkCellReadsStats.awk CellReads.stats Features.stats -v flags=mult

BEGIN {
    OFS="\t";
    if (flags~"mult") {
        flagMult=1;
        print "Yes multi-gene";
    } else {
        print "No multi-gene";
    };
}

(ARGIND==1) {
if (NR==1) {
    for (ff=1; ff<=NF; ff++) {
        fn[ff]=$ff;
        fi[$ff]=ff;
    };
    next;
};

if ($1~/^CB/) {
    for (ff=2; ff<=NF; ff++) {
        nocbStats[ff] += $ff;
    };
    next;
};
# check fields
if ($fi["cbMatch"] != $fi["cbPerfect"]+$fi["cbMMunique"]+$fi["cbMMmultiple"]) {
    print;
    print "cbMatch != cbPerfect + cbMMunique + cbMMmultiple";
    exit;
};
if ( flagMult && $fi["featureM"]!=$fi["countedM"] ) {
    print;
    print "featureM!=countedM";
    exit;
};
if ( $fi["featureU"] != $fi["countedU"] ) {
    print;
    print "featureU!=countedU";
    exit;
};

for (ff=2; ff<=NF; ff++) {
    sumStats[ff] += $ff;
};

}

(ARGIND==2) {
    featStats[$1]=$2;
}

END {
   #print "featureM", sumStats[fi["featureM"]], featStats["MultiFeature"];
   # featureM counts reads with defined CB, while MultiFeature also counts reads that may be CB-rejected in inputRecords
   print "countedU+M", sumStats[fi["countedU"]]+sumStats[fi["countedM"]], featStats["yesWLmatch"];
   print "countedU", sumStats[fi["countedU"]], featStats["yessubWLmatch_UniqueFeature"];
   print "nUMIunique", sumStats[fi["nUMIunique"]], featStats["yesUMIs"];

   # this does not work: cbMatch (all cells with valid CBs) cannot be extracted from Features.stats, 
   # because some of the cells with possible valid barcodes will be classified as noFeature,unmapped,multifeature
   #split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches MultiFeature yesWLmatch", tags1);
   #split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches yesWLmatch", tags1);
   #split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches MultiFeature yessubWLmatch_UniqueFeature", tags1);
   #for (tt in tags1) {
   #    sum1+=featStats[tags1[tt]];
   #};
   #print "cbMatch", sumStats[fi["cbMatch"]], sum1;

   print "---------------------------------- All Sums";
   print "allReads", nocbStats[fi["cbMatch"]]+sumStats[fi["cbMatch"]];
   for (f=2;f<=length(fn);f++)
       print fn[f], sumStats[f]+0;
}
