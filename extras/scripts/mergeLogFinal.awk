#
# merges Log.final.out files from multiple runs into one table
# usage:
# awk -f mergeLogFinal.awk /path/to/1st/Log.final.out /path/to/2nd/Log.final.out ...
# e.g.
# awk -f mergeLogFinal.awk */Log.final.out
#

BEGIN {
    FS="|";
    for (jj=1;jj<=ARGC;jj++)
    {
        a=ARGV[jj]; 
        gsub("/Log.final.out","",a);
        printf ";" a
    }; 
    printf "\n";
} 
{
    gsub(/^[ \t]+|[ \t]+$/,"",$1);
    gsub(/^[ \t]+|[ \t]+$/,"",$2);
    L[FNR]=$1; 
    V[FNR,ARGIND]=$2
} 
END {
    for (ii=1;ii<=length(L);ii++) 
    {
        printf "%s",L[ii];
        if (V[ii,1]!="") 
            for (jj=1;jj<=ARGC;jj++) 
                printf ";" V[ii,jj]; 
        printf "\n"
     } 
}
