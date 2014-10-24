BEGIN {
  nSection=0;
  optOptTableEnd="\\end{optOptTable}";
  optOptTableBegin="\\begin{optOptTable}";
  optTableEnd="\\end{optTable}";
  optTableBegin="\\begin{optTable}";
};

{
if ($1=="###") {# new group/subsection of parameters
  if ($2!="versions") {# skip versions
      if (nSection>0) print optTableEnd;
      print "\\optSection{" substr($0,length($1)+2,1) tolower(substr($0,length($1)+3)) "}";
      print optTableBegin;
      ++nSection;
  };
} else if ($0!="" && substr($0,1,1)!=" " && substr($1,1,1)!="#" && substr($1,1,7)!="version") {//option name has a letter as the first character
  optV=$2; 
  for (ii=3;ii<=NF;ii++) optV=" " $ii;
  print "\\optName{" $1 "}" " & " "\\optValue{" optV "}" " \\\\";

  getline;
  nOpt=0;
  while ($0!="") {
      $0=substr($0,match($0,/[^[:space:]]/));
      no=split($0,oo,/[:space:]*\.\.\.[:space:]*/);
      if (no!=2) {# not option line
           if (nOpt>0) print optOptTableEnd "\\\\";
           print " & " "\\optLine{" $0 "}" "\\\\";
           nOpt=0;
      } else {
           if (nOpt==0) print " & " optOptTableBegin;
           #n1=match($0,/[:space:]*\.\.\./);
           print "\\optOpt{" substr(oo[1],1,match(oo[1],/[[:space:]$]/)-1) "}" " & " "\\optOptLine{" substr(oo[2],match(oo[2],/[^[:space:]]/)) "}" "\\\\";
           ++nOpt;
      };
      getline;
  };
 if (nOpt>0) print optOptTableEnd "\\\\";
};

};

END {
    print optTableEnd;
};
