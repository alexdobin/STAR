###
# usage
# awk -f mergeSuperContig.awk All.fasta All.gtf
###

BEGIN {
    ################################
    # these parameters can be changed
    
    shortL=64000
    
    fLongOut="Long.out.fasta";     # output file name for FASTA
    fShortOut="Short.out.fasta";   # output file name for FASTA
    fGtfOut="Annot.out.gtf";       # output file name for GTF
    fChrStartOut="ChrStart.tab";   # output file name for chr starts
    pN=60;                         # N-padding length
    scName="SC";                   # supercontig name
    ################################

    OFS="\t";
    s=0;  # sequence start in super-contig
    c=""; # current chr name
    Seq=""; # chr sequence

    for (ii=1; ii<=60; ii++) P=P "N"; # N-padding string

    print ">" scName > fShortOut;

    while (getline < ARGV[1])
    {# FASTA
        if (substr($1,1,1)==">")
        {# new chr
           chrFinish();
           # start new chr
           c=substr($1,2); # chr name
           l=0;            # chr length
           Seq="";
           print c,s > fChrStartOut;
        } else
        {# collect sequence
           Seq=Seq $1;
           l+=length($1);
        };          
     };

     chrFinish();

     while (getline < ARGV[2])
     {# GTF
          if ($1 in S)
          {# transform coordinates
              $4+=S[$1];
              $5+=S[$1];
              $1=scName;
          };
      print > fGtfOut;  
      };         
};

function chrFinish()
{
           if (l>shortL)
           {# long chr
               print ">" c "\n" Seq > fLongOut; 
           } else if (l>0)
           {#short chr
               printf Seq P > fShortOut;
               S[c]=s;         # chr start in contig
               s+=l+pN;
           };
};
