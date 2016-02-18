#include "Stats.h"
#include "TimeFunctions.h"

void Stats::resetN() {//zero all counters
    readN = 0; readBases = 0;
    mappedMismatchesN = 0; mappedInsN = 0; mappedDelN = 0; mappedInsL = 0; mappedDelL = 0; mappedBases = 0;  mappedPortion = 0;
    mappedReadsU = 0; mappedReadsM = 0;
    unmappedOther = 0; unmappedShort = 0; unmappedMismatch = 0; unmappedMulti = 0; unmappedAll = 0;
    chimericAll = 0;
    splicesNsjdb=0;
    for (uint ii=0; ii<SJ_MOTIF_SIZE; ii++) {
        splicesN[ii]=0;
    };

};

Stats::Stats() {//constructor
    resetN();
    timeLastReport=0;
};

void Stats::addStats(Stats &S) {//add S to Stats
    readN += S.readN; readBases += S.readBases;
    mappedMismatchesN += S.mappedMismatchesN; mappedInsN += S.mappedInsN; mappedDelN += S.mappedDelN;
    mappedInsL += S.mappedInsL; mappedDelL += S.mappedDelL; mappedBases += S.mappedBases;  mappedPortion += S.mappedPortion;
    mappedReadsU += S.mappedReadsU; mappedReadsM += S.mappedReadsM;
    unmappedOther += S.unmappedOther; unmappedShort += S.unmappedShort; unmappedMismatch += S.unmappedMismatch; unmappedMulti += S.unmappedMulti; unmappedAll += S.unmappedAll;
    chimericAll += S.chimericAll;

    splicesNsjdb += S.splicesNsjdb;
    for (uint ii=0; ii<SJ_MOTIF_SIZE; ii++) {
        splicesN[ii] +=S.splicesN[ii];
    };
};

void Stats::transcriptStats(Transcript &T, uint Lread) {
    mappedMismatchesN += T.nMM;
    mappedInsN += T.nIns;
    mappedDelN += T.nDel;
    mappedInsL += T.lIns;
    mappedDelL += T.lDel;

    uint mappedL=0;
    for (uint ii=0; ii<T.nExons; ii++) {
        mappedL += T.exons[ii][EX_L];
    };
    for (uint ii=0; ii<T.nExons-1; ii++) {
        if (T.canonSJ[ii]>=0) splicesN[T.canonSJ[ii]]++;
        if (T.sjAnnot[ii]==1) splicesNsjdb++;
    };

    mappedBases += mappedL;
    mappedPortion += double(mappedL)/double(Lread);
};

#define SETW1 setw(9)
#define SETW2 setw(8)
#define SETW3 setw(12)

void Stats::progressReportHeader(ofstream &progressStream) {
    progressStream  <<setw(15)<< "Time" <<SETW1<< "Speed" <<SETW3<< "Read" <<SETW1<< "Read" <<SETW1<< "Mapped" \
            <<SETW1<< "Mapped" <<SETW1<< "Mapped" \
            <<SETW1<< "Mapped" <<SETW1<< "Unmapped" <<SETW1<< "Unmapped" <<SETW1<< "Unmapped" <<SETW1<< "Unmapped"\
            << "\n";
    progressStream  <<setw(15)<< " " <<SETW1<< "M/hr" <<SETW3<< "number" <<SETW1<< "length" <<SETW1<< "unique" \
            <<SETW1<< "length" <<SETW1<< "MMrate" \
            <<SETW1<< "multi" <<SETW1<< "multi+" <<SETW1<< "MM" <<SETW1<< "short" <<SETW1<< "other"\
            << "\n"<<flush;
};

void Stats::progressReport(ofstream &progressStream) {

    time_t timeCurrent;
    time( &timeCurrent);

    if (difftime(timeCurrent,timeLastReport)>=60.0 && readN>0) {//make the report
        //progressStream.imbue(std::locale(""));
        progressStream <<setw(15)<< timeMonthDayTime(timeCurrent) \
                <<SETW1<< setiosflags(ios::fixed) << setprecision(1) \
                << double(readN)/1e6/difftime(timeCurrent,timeStartMap)*3600 \
                <<SETW3<< readN \
                <<SETW1<< (readN>0 ? readBases/readN : 0) \
                <<SETW2<< (readN>0 ? double(mappedReadsU)/double(readN)*100 : 0)  <<'%' \
                <<SETW1<< (readN>0 ? double(mappedBases)/double(mappedReadsU) : 0)
                <<SETW2<< (readN>0 ? double(mappedMismatchesN)/double(mappedBases)*100 : 0) <<'%' \
                <<SETW2<< (readN>0 ? double(mappedReadsM)/double(readN)*100 : 0) <<'%'\
                <<SETW2<< (readN>0 ? double(unmappedMulti)/double(readN)*100 : 0) <<'%'\
                <<SETW2<< (readN>0 ? double(unmappedMismatch)/double(readN)*100 : 0) <<'%'\
                <<SETW2<< (readN>0 ? double(unmappedShort)/double(readN)*100 : 0)<<'%'\
                <<SETW2<< (readN>0 ? double(unmappedOther)/double(readN)*100 : 0) <<'%'\
                <<"\n"<<flush;
        timeLastReport=timeCurrent;

    };
};

void Stats::reportFinal(ofstream &streamOut) {
    int w1=50;
    time( &timeFinish);

                                            //<<setiosflags(ios::left)
    streamOut  <<setiosflags(ios::fixed)  << setprecision(2) \
               <<setw(w1)<< "Started job on |\t" << timeMonthDayTime(timeStart)<<"\n" \
               <<setw(w1)<< "Started mapping on |\t" << timeMonthDayTime(timeStartMap)<<"\n" \
               <<setw(w1)<< "Finished on |\t"<< timeMonthDayTime(timeFinish)<<"\n" \
               <<setw(w1)<< "Mapping speed, Million of reads per hour |\t"<< double(readN)/1e6/difftime(timeFinish,timeStartMap)*3600<<"\n" \
               <<"\n" \
               <<setw(w1)<< "Number of input reads |\t"                        << readN <<"\n" \
               <<setw(w1)<< "Average input read length |\t"                    << (readN>0 ? readBases/readN : 0) <<"\n" \
               <<setw(w1)<< "UNIQUE READS:\n" \
               <<setw(w1)<< "Uniquely mapped reads number |\t"                 << mappedReadsU <<"\n" \
               <<setw(w1)<< "Uniquely mapped reads % |\t"                      << (readN>0 ? double(mappedReadsU)/double(readN)*100 : 0) <<'%'<<"\n" \
               <<setw(w1)<< "Average mapped length |\t"                        << (mappedReadsU>0 ? double(mappedBases)/double(mappedReadsU) : 0) <<"\n";

    streamOut  <<setw(w1)<< "Number of splices: Total |\t"                     << splicesN[0]+splicesN[1]+splicesN[2]+splicesN[3]+splicesN[4]+splicesN[5]+splicesN[6]<< "\n" \
               <<setw(w1)<< "Number of splices: Annotated (sjdb) |\t"          << splicesNsjdb << "\n" \
               <<setw(w1)<< "Number of splices: GT/AG |\t"                     << splicesN[1]+splicesN[2] << "\n" \
               <<setw(w1)<< "Number of splices: GC/AG |\t"                     << splicesN[3]+splicesN[4] << "\n" \
               <<setw(w1)<< "Number of splices: AT/AC |\t"                     << splicesN[5]+splicesN[6] << "\n" \
               <<setw(w1)<< "Number of splices: Non-canonical |\t"             << splicesN[0] << "\n";

    streamOut  <<setw(w1)<< "Mismatch rate per base, % |\t"                << double(mappedMismatchesN)/double(mappedBases)*100 <<'%' <<"\n" \
               <<setw(w1)<< "Deletion rate per base |\t"                       << (mappedBases>0 ? double(mappedDelL)/double(mappedBases)*100 : 0) <<'%' <<"\n" \
               <<setw(w1)<< "Deletion average length |\t"                      << (mappedDelN>0 ? double(mappedDelL)/double(mappedDelN) : 0) <<"\n" \
               <<setw(w1)<< "Insertion rate per base |\t"                      << (mappedBases>0 ? double(mappedInsL)/double(mappedBases)*100 : 0) <<'%' <<"\n" \
               <<setw(w1)<< "Insertion average length |\t"                     << (mappedInsN>0 ? double(mappedInsL)/double(mappedInsN) : 0) <<"\n" \
               <<setw(w1)<< "MULTI-MAPPING READS:\n" \
               <<setw(w1)<< "Number of reads mapped to multiple loci |\t"      << mappedReadsM <<"\n" \
               <<setw(w1)<< "% of reads mapped to multiple loci |\t"           << (readN>0 ? double(mappedReadsM)/double(readN)*100 : 0)<<'%' <<"\n" \
               <<setw(w1)<< "Number of reads mapped to too many loci |\t"      << unmappedMulti <<"\n" \
               <<setw(w1)<< "% of reads mapped to too many loci |\t"           << (readN>0 ? double(unmappedMulti)/double(readN)*100 : 0) <<'%' <<"\n" \
               <<setw(w1)<< "UNMAPPED READS:\n" \
               <<setw(w1)<< "% of reads unmapped: too many mismatches |\t"     << (readN>0 ? double(unmappedMismatch)/double(readN)*100 : 0) <<'%' <<"\n" \
               <<setw(w1)<< "% of reads unmapped: too short |\t"               << (readN>0 ? double(unmappedShort)/double(readN)*100 : 0) <<'%' <<"\n" \
               <<setw(w1)<< "% of reads unmapped: other |\t"                   << (readN>0 ? double(unmappedOther)/double(readN)*100 :0) <<'%'<<"\n" \
               <<setw(w1)<< "CHIMERIC READS:\n" \
               <<setw(w1)<< "Number of chimeric reads |\t"                     << chimericAll <<"\n" \
               <<setw(w1)<< "% of chimeric reads |\t"                          << (readN>0 ? double(chimericAll)/double(readN)*100 :0) <<'%'<<"\n" <<flush;

};


