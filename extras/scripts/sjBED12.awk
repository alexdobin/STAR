# for SJs (e.g.converts 3-column BED into BED12
# awk -v blockLen=10 -f ... sjIn.bed 

BEGIN {
    OFS="\t";

}

{
    if ($1=="")
        next;
    #################### input
    chrom=$1;
    chromStart1=$2;
    chromEnd1=$3;
    name="sj-" NR;
    if (name=="")
        name=".";

    score=1000;

    strand = $6;
    if (strand!="+" && strand!="-")
        strand=".";

    itemRgb="180,0,0";

    #################### BED12
    chromStart = chromStart1-blockLen;
    chromEnd = chromEnd1+blockLen;
    thickStart = chromStart;
    thickEnd = chromEnd;
    blockCount = 2;
    blockSizes = blockLen "," blockLen;
    blockStarts = 0 "," chromEnd1-chromStart;

    print chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts > FILENAME ".bed12";

    ####### GTF
    source = "STAR_SJ";
    feature = "exon";
    frame = ".";
    group = "transcript_id \"" name "\"; gene_id \"" name "\";" 

    print chrom, source, feature, chromStart1-blockLen+1, chromStart1, score, strand, frame, group > FILENAME ".gtf";
    print chrom, source, feature, chromEnd1+1, chromEnd1+blockLen, score, strand, frame, group > FILENAME ".gtf";
}
