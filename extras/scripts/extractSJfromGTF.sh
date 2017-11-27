# usage:
# extractSJfromGTF.sh in.gtf > out.sj
#
# assumes transcript_id in the 12th field of GTF
#
awk '$3=="exon" {print $12,$1,$4,$5,$7}' $1 |\
 sort -k1,1V -k2,2V -k3,3n |\
 awk 'BEGIN {OFS="\t"} {if (t==$1) {print $2,e1+1,$3-1,$5}; e1=$4;t=$1 }' |\
 sort -k1,1V -k2,2n -k3,3n -k4,4 | uniq

