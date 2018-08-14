#
# merges ReadsPerGene.out.tab files from multiple runs into one table
# usage:
# awk -f mergeGeneCounts.awk -v Col=2 /path/to/1st/ReadsPerGene.out.tab /path/to/2nd/ReadsPerGene.out.tab ...
# e.g.
# awk -f mergeGeneCounts.awk -v Col=2 */ReadsPerGene.out.tab
#
# -v Col=<column to add to the table>: depends on the standedness of the table
# advanced parameters
# -v Skip=<number of lines to skip>
# -v Name=<common file name substring to remove from column names>




BEGIN {
    FS="\t";

    if (Name=="") Name="/ReadsPerGene.out.tab";
    if (Skip=="") Skip=0;
    if (Col=="") {
        print "Specify the column with -v Col=..." > /dev/err;
        exit;
    };

    for (jj=1;jj<=ARGC;jj++)
    {# print header line with file names
        a=ARGV[jj]; 
        gsub(Name,"",a);
        printf ";" a
    }; 
    printf "\n";
} 
{
    if (ARGIND==1) {
        L[FNR]=$1; # record gene names (1st column)
    } else {
        if ($1!=L[FNR]) {
            print "File #" ARGIND ": " FILENAME " is not sorted properly, sort all files by the first column" >/dev/err;
        };
    };

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
