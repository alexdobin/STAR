# usage:
# cat Aligned.out.sam | awk -v strType=2 -f tagXSstrandedData.awk 
# strType defines strandedness of the libraries: strType = mate whose strand is the same as RNA strand.
# For instance, for Illumina Tru-seq, strType=2 - the 2nd mate's strand is the same as RNA.


BEGIN {
    OFS="\t";
    strSym[0]="+";
    strSym[1]="-";
}

{

    if (substr($1,1,1)=="@" || $4==0)
    {# header, or unmapped read - just print
        print;
        next;
    };

    str=and($2,0x10)/0x10;
   
    if (and($2,0x1)==0)
    {# single end defaults to mate
        mate=1;
    } else
    {
        mate=and($2,0x40)/0x40+2*and($2,0x80)/0x80;
    };

    if (mate>0 && mate <3)
    {# mate is defined - add XS tag
       if (mate!=strType) str=1-str; #revert strand if the mate is opposite
       print $0 "\t" "XS:A:" strSym[str];
    } else 
    {# mate is not defined - just print
       print;
    };    
}
