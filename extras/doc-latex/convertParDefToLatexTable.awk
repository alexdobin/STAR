function substLatexSymbols() {
 gsub("\\^","\\^{}",$0);
 gsub(">","{\\textgreater}");
 gsub("<","{\\textless}");
 gsub("_","{\\textunderscore}");
 gsub("&","{\\\\&}");
};

BEGIN {
  nSection=0;
  optOptTableEnd="\\end{optOptTable}";
  optOptTableBegin="\\begin{optOptTable}";
  optTableEnd="\\end{optTable}";
  optTableBegin="\\begin{optTable}";
};

{
substLatexSymbols();
if ($1=="###") {# new group/subsection of parameters
  if ($2!="versions") {# skip versions
      if (nSection>0) print optTableEnd;
      sectionName=substr($0,index($0,$2));
      printf "\\optSection{" sectionName "}";
      gsub(" ","_",sectionName);
      print  "\\label{" sectionName "}";
      print optTableBegin;
      ++nSection;
  };
} else if ($0!="" && substr($0,1,1)!=" " && substr($1,1,1)!="#" && substr($1,1,7)!="version") {//option name has a letter as the first character
  optV=$2; 
  for (ii=3;ii<=NF;ii++) optV=optV " " $ii;
  print "\\optName{" $1 "}";
  print "  \\optValue{" optV "}";

  getline;substLatexSymbols();
  nOpt=0;
  while ($1!="") {
      $0=substr($0,match($0,/[^[:space:]]/));
      no=split($0,oo,/[:space:]*\.\.\.[:space:]*/);
      if (no!=2) {# not option line
           if (nOpt>0) print optOptTableEnd;
           print "  \\optLine{" $0 "}" " ";
           nOpt=0;
      } else {
           if (nOpt==0) print optOptTableBegin;
           gsub(/^[ \t]+|[ \t]+$/, "",oo[1]); #remove leading trailing spaces
           gsub(/^[ \t]+|[ \t]+$/, "",oo[2]);
           print "  \\optOpt{" oo[1] "}   \\optOptLine{" oo[2] "}" ;
           ++nOpt;
      };
      getline; substLatexSymbols();
  };
  if (nOpt>0) print optOptTableEnd;
};

};

END {
    print optTableEnd;
};
