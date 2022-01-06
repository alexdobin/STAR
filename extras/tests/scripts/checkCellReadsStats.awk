# usage awk -f checkCellReadsStats.awk CellReads.stats Features.stats

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

if ($1~/^CB/)
    next;
# check fields
if ($fi["cbMatch"] != $fi["cbPerfect"]+$fi["cbMMunique"]+$fi["cbMMmultiple"]) {
    print;
    print "cbMatch != cbPerfect + cbMMunique + cbMMmultiple";
    exit;
};
if ($fi["featureU"] != $fi["countedU"] || $fi["featureM"]!=$fi["countedM"]) {
    print;
    print "feature!=counted";
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
   print "featureM", sumStats[fi["featureM"]], featStats["MultiFeature"];
   print "countedU+M", sumStats[fi["countedU"]]+sumStats[fi["countedM"]], featStats["yesWLmatch"];
   print "countedU", sumStats[fi["countedU"]], featStats["yessubWLmatch_UniqueFeature"];
   print "nUMIunique", sumStats[fi["nUMIunique"]], featStats["yesUMIs"];

   #split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches MultiFeature yesWLmatch", tags1);
   split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches yesWLmatch", tags1);
   #split("noUnmapped noNoFeature noMMtoWLwithoutExact noTooManyWLmatches MultiFeature yessubWLmatch_UniqueFeature", tags1);
   for (tt in tags1) {
       sum1+=featStats[tags1[tt]];
   };
   print "cbMatch", sumStats[fi["cbMatch"]], sum1;

   print "---------------------------------- All Sums";
   for (f=2;f<=length(fn);f++)
       print fn[f], sumStats[f]+0;
}
