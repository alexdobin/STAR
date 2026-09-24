#!/bin/bash
# Regression test for transcripts/genes whose GTF strand is undefined ('.'):
#  1. Aligned.toTranscriptome.out.bam records of a '.' transcript must have the same flags and
#     sequences as those of an identical '+' transcript (before the fix the strand flag was
#     inverted and the sequence reverse-complemented at unchanged + coordinates);
#  2. STARsolo --soloStrand Forward must count reads of both strands for a '.' gene (before the
#     fix only antisense reads were counted).
# Usage: testUndefinedStrand.sh [/path/to/STAR]   (needs samtools in PATH for check 1)
# Exit status 0 = pass.
set -e
STAR=${1:-$(dirname "$0")/../../../source/STAR}
T=$(mktemp -d)
cd "$T"
# 3-kb pseudo-random genome (LCG in awk, implementation-independent), two identical two-exon
# gene structures: TP on '+' (exons 201-500, 801-1100) and TD on '.' (1701-2000, 2301-2600),
# GT..AG planted at both introns.
awk 'BEGIN{x=12345; s=""; for(i=0;i<3000;i++){x=(x*1103515245+12345)%2147483648; s=s substr("ACGT",int(x/65536)%4+1,1)};
  s=substr(s,1,499) "CGT" substr(s,503,296) "AGC" substr(s,802);
  s=substr(s,1,1999) "CGT" substr(s,2003,296) "AGC" substr(s,2302);
  print ">chr1"; print s}' > g.fa
printf 'chr1\tt\texon\t201\t500\t.\t+\t.\tgene_id "gTP"; transcript_id "TP";\nchr1\tt\texon\t801\t1100\t.\t+\t.\tgene_id "gTP"; transcript_id "TP";\n' > g.gtf
printf 'chr1\tt\texon\t1701\t2000\t.\t.\t.\tgene_id "gTD"; transcript_id "TD";\nchr1\tt\texon\t2301\t2600\t.\t.\t.\tgene_id "gTD"; transcript_id "TD";\n' >> g.gtf
# reads: exonic (offset 100) and junction-spanning (offset 270) 60-mers of each transcript, forward and reverse-complemented
awk 'NR==2{g=$0; n=split("TP 201 500 801 1100,TD 1701 2000 2301 2600",tr,",");
  for(i=1;i<=n;i++){split(tr[i],f," "); t=substr(g,f[2],f[3]-f[2]+1) substr(g,f[4],f[5]-f[4]+1);
    for(o=100;o<=270;o+=170){r=substr(t,o+1,60); rc=""; for(k=60;k>=1;k--){c=substr(r,k,1); rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")};
      print "@" f[1] "_fwd_" o "\n" r "\n+\n" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
      print "@" f[1] "_rev_" o "\n" rc "\n+\n" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}}}' g.fa > r.fq
# barcode reads for STARsolo: one cell barcode, non-homopolymer UMIs at least 6 mismatches apart (not collapsed by --soloUMIdedup 1MM_All, not removed as homopolymers)
awk 'BEGIN{split("ACACACACACAC AGAGAGAGAGAG ATATATATATAT CGCGCGCGCGCG CTCTCTCTCTCT GTGTGTGTGTGT AACCAACCAACC GGTTGGTTGGTT",u," ")} NR%4==1{print $0; print "ACGTACGTACGTACGT" u[(NR-1)/4+1]; print "+"; print "IIIIIIIIIIIIIIIIIIIIIIIIIIII"}' r.fq > b.fq
echo ACGTACGTACGTACGT > wl.txt
mkdir gidx
"$STAR" --runMode genomeGenerate --genomeDir gidx --genomeFastaFiles g.fa --sjdbGTFfile g.gtf --sjdbOverhang 59 --genomeSAindexNbases 4 --outFileNamePrefix gidx/ > /dev/null
"$STAR" --genomeDir gidx --readFilesIn r.fq --outFileNamePrefix tr/ --outSAMtype None --quantMode TranscriptomeSAM > /dev/null
"$STAR" --genomeDir gidx --readFilesIn r.fq b.fq --outFileNamePrefix solo/ --outSAMtype None --soloType CB_UMI_Simple --soloCBwhitelist wl.txt --soloUMIlen 12 \
        --soloFeatures Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS --soloStrand Forward > /dev/null
fail=0
if command -v samtools > /dev/null; then
    # flag, position and CIGAR of each TD record must equal those of the TP record with the same suffix
    samtools view tr/Aligned.toTranscriptome.out.bam | awk '{sub(/^T[PD]_/,"",$1); print $1, $2, $4, $6}' | sort | uniq -c | awk '$1!=2{bad=1; print "MISMATCH between TP and TD records (name flag pos CIGAR):", $0} END{exit bad}' || fail=1
else
    echo "samtools not found: transcriptome BAM flags not checked"
fi
for f in Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS; do
    n=$(awk 'NR>3 && $1==2{s+=$3} END{print s+0}' solo/Solo.out/$f/raw/matrix.mtx)
    m=$(awk 'NR>3 && $1==1{s+=$3} END{print s+0}' solo/Solo.out/$f/raw/matrix.mtx)
    if [ "$n" != 4 ] || [ "$m" != 2 ]; then echo "STARsolo $f --soloStrand Forward: gTD counted $n of 4 reads (expected 4), gTP $m (expected 2)"; fail=1; fi
done
cd /; rm -rf "$T"
if [ $fail = 0 ]; then echo "testUndefinedStrand: PASS"; else echo "testUndefinedStrand: FAIL"; fi
exit $fail
