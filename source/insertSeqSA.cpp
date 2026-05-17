/*
 * inserts sequences into the SA
 * returns number of SA indexes inserted
 */
#include "insertSeqSA.h"
#include "ErrorWarning.h"
#include "SuffixArrayFuns.h"
#include "GenomeInsertOverlay.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include "streamFuns.h"
#include "binarySearch2.h"
#include "funCompareUintAndSuffixes.h"
#include "funCompareUintAndSuffixesMemcmp.h"
#include <cmath>
#include <map>
#include "sortSuffixesBucket.h"

namespace {
const string deltaMagic = "STARgenomeInsertDeltaV1";

uint64 genomeInsertDeltaHash(const char *seq, uint64 length)
{
    uint64 hash=14695981039346656037ULL;
    for (uint64 ii=0; ii<length; ii++) {
        hash ^= (uint64) (unsigned char) seq[ii];
        hash *= 1099511628211ULL;
    };
    return hash;
}

uint64 genomeInsertCountBefore(uint64 *indArray, uint64 nInd, uint64 oldSAindex)
{
    uint64 left=0;
    uint64 right=nInd;
    while (left<right) {
        uint64 middle=(left+right)/2;
        if (indArray[2*middle]<oldSAindex) {
            left=middle+1;
        } else {
            right=middle;
        };
    };
    return left;
}

uint genomeInsertInsertedSAvalue(uint64 insertIndex, uint64 nG, uint64 nG1, uint64 nG2, uint N2bit)
{
    if (insertIndex<nG1) {
        return insertIndex+nG;
    } else {
        return (insertIndex-nG1+nG2) | N2bit;
    };
}

bool genomeInsertOutputByteAligned(uint64 outputIndex, uint wordLength)
{
    return ((outputIndex*wordLength) % 8)==0;
}

uint64 genomeInsertOutputIndex(uint64 oldSAindex, uint64 *indArray, uint64 nInd)
{
    return oldSAindex + genomeInsertCountBefore(indArray, nInd, oldSAindex);
}

uint64 genomeInsertFindAlignedBoundary(uint64 oldSAindex, uint64 oldSAend, uint64 *indArray, uint64 nInd, uint wordLength)
{
    for (uint64 old1=oldSAindex; old1<oldSAend; old1++) {
        if (genomeInsertOutputByteAligned(genomeInsertOutputIndex(old1, indArray, nInd), wordLength)) {
            return old1;
        };
    };
    return oldSAend;
}

void writePackedToZeroedBytes(char *charArray, uint wordLength, uint64 index, uint value)
{
    // Chunk writers own disjoint byte-aligned ranges, so byte-wise OR avoids
    // the overlapping unaligned writes used by PackedArray::writePacked.
    unsigned char *byteArray=(unsigned char*) charArray;
    uint64 bitIndex=index*wordLength;
    uint64 byteIndex=bitIndex/8;
    uint bitShift=bitIndex%8;
    uint valueMask = wordLength>=sizeof(uint)*8 ? (uint) -1 : ((1LLU << wordLength) - 1);
    uint valueShifted=(value & valueMask) << bitShift;
    uint nBytes=(wordLength+bitShift+7)/8;

    for (uint ii=0; ii<nBytes; ii++) {
        byteArray[byteIndex+ii] |= (unsigned char) ((valueShifted >> (8*ii)) & 0xFF);
    };
}

void insertSeqSAwriteChunk(PackedArray &SA1, vector<uint> &oldSAvalues, uint64 oldValuesBegin, uint wordLength, uint64 *indArray, uint64 nInd, uint64 oldBegin, uint64 oldEnd, uint64 nG, uint64 nG1, uint64 nG2, uint N2bit)
{
    uint64 insertBegin=genomeInsertCountBefore(indArray, nInd, oldBegin);
    uint64 insertEnd=genomeInsertCountBefore(indArray, nInd, oldEnd);
    uint64 outputBegin=oldBegin+insertBegin;
    uint64 outputEnd=oldEnd+insertEnd;
    uint64 byteBegin=outputBegin*wordLength/8;
    uint64 byteEnd=outputEnd*wordLength/8;
    char *chunkOut=SA1.charArray+byteBegin;

    memset(chunkOut, 0, byteEnd-byteBegin);

    uint64 insert1=insertBegin;
    uint64 output1=0;
    for (uint64 old1=oldBegin; old1<oldEnd; old1++) {
        while (insert1<insertEnd && old1==indArray[2*insert1]) {
            writePackedToZeroedBytes(chunkOut, wordLength, output1, genomeInsertInsertedSAvalue(indArray[2*insert1+1], nG, nG1, nG2, N2bit));
            ++output1;
            ++insert1;
        };
        writePackedToZeroedBytes(chunkOut, wordLength, output1, genomeInsertSAshift(oldSAvalues[old1-oldValuesBegin], nG, nG1, nG2, wordLength-1));
        ++output1;
    };
}

uint64 genomeInsertDeltaHeaderValue(const map<string,uint64> &header, const string &key, const string &fileName, Parameters &P)
{
    map<string,uint64>::const_iterator value=header.find(key);
    if (value==header.end()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: genome insert delta file is missing header value " << key << "\n";
        errOut << "Delta file: " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    return value->second;
}

void writeGenomeInsertDelta(const string &fileName, uint64 *indArray, uint64 nInd, uint64 nG, uint64 nG1, uint64 nG2, const char *G1, PackedArray &SA, Parameters &P, Genome &mapGen)
{
    ofstream deltaOut(fileName.c_str(), ios::out | ios::binary);
    if (deltaOut.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not create genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
    };

    deltaOut << deltaMagic << "\n";
    deltaOut << "nG\t" << nG << "\n";
    deltaOut << "nG1\t" << nG1 << "\n";
    deltaOut << "nG2\t" << nG2 << "\n";
    deltaOut << "nSA\t" << SA.length << "\n";
    deltaOut << "SAwordLength\t" << SA.wordLength << "\n";
    deltaOut << "gSAindexNbases\t" << mapGen.pGe.gSAindexNbases << "\n";
    deltaOut << "gSuffixLengthMax\t" << mapGen.pGe.gSuffixLengthMax << "\n";
    deltaOut << "insertSeqHash\t" << genomeInsertDeltaHash(G1, nG1) << "\n";
    deltaOut << "nInd\t" << nInd << "\n";
    deltaOut << "END\n";
    deltaOut.write((char*) indArray, (streamsize) (2*nInd*sizeof(uint64)));
    if (deltaOut.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error while writing genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_FILE_WRITE, P);
    };
    deltaOut.close();
};

uint64 readGenomeInsertDelta(const string &fileName, uint64 *indArray, uint64 indArrayCapacity, uint64 nG, uint64 nG1, uint64 nG2, const char *G1, PackedArray &SA, Parameters &P, Genome &mapGen)
{
    ifstream deltaIn(fileName.c_str(), ios::in | ios::binary);
    if (deltaIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: could not open genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    string line;
    getline(deltaIn, line);
    if (line!=deltaMagic) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: unsupported genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    map<string,uint64> header;
    while (getline(deltaIn, line)) {
        if (line=="END") {
            break;
        };
        istringstream lineStream(line);
        string key;
        uint64 value=0;
        if (!(lineStream >> key >> value)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal INPUT FILE error: malformed genome insert delta header in " << fileName << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        header[key]=value;
    };
    if (line!="END") {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: malformed genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (genomeInsertDeltaHeaderValue(header, "nG", fileName, P)!=nG
            || genomeInsertDeltaHeaderValue(header, "nG1", fileName, P)!=nG1
            || genomeInsertDeltaHeaderValue(header, "nG2", fileName, P)!=nG2
            || genomeInsertDeltaHeaderValue(header, "nSA", fileName, P)!=SA.length
            || genomeInsertDeltaHeaderValue(header, "SAwordLength", fileName, P)!=SA.wordLength
            || genomeInsertDeltaHeaderValue(header, "gSAindexNbases", fileName, P)!=mapGen.pGe.gSAindexNbases
            || genomeInsertDeltaHeaderValue(header, "gSuffixLengthMax", fileName, P)!=mapGen.pGe.gSuffixLengthMax
            || genomeInsertDeltaHeaderValue(header, "insertSeqHash", fileName, P)!=genomeInsertDeltaHash(G1, nG1)) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: genome insert delta file does not match the loaded base genome/index or inserted FASTA files\n";
        errOut << "Delta file: " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    uint64 nInd=genomeInsertDeltaHeaderValue(header, "nInd", fileName, P);
    if (nInd>indArrayCapacity) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error: genome insert delta file contains too many suffix indices: " << nInd << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    deltaIn.read((char*) indArray, (streamsize) (2*nInd*sizeof(uint64)));
    if (deltaIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal INPUT FILE error while reading genome insert delta file " << fileName << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    return nInd;
};

// Update the prefix table for newly inserted suffixes without rebuilding it from the expanded SA.
void insertSeqSAi(PackedArray &SAi, char **seq1, uint64 *indArray, uint64 nInd, Parameters &P, Genome &mapGen)
{
    time_t rawtime;

    for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {
        uint64 iSeq=0;
        uint ind0=mapGen.genomeSAindexStart[iL]-1;//last index that was present in the old genome

        for (uint ii=mapGen.genomeSAindexStart[iL]; ii<mapGen.genomeSAindexStart[iL+1]; ii++) {
            uint iSA1=SAi[ii];
            uint iSA2=iSA1 & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask;

            if (iSeq<nInd && (iSA1 & mapGen.SAiMarkAbsentMaskC)>0) {
                uint64 iSeq1=iSeq;
                int64 ind1=funCalcSAi(seq1[0]+indArray[2*iSeq+1], iL);
                while (iSeq<nInd && ind1 < (int64)(ii-mapGen.genomeSAindexStart[iL]) && indArray[2*iSeq]<=iSA2) {
                    ++iSeq;
                    if (iSeq<nInd) {
                        ind1=funCalcSAi(seq1[0]+indArray[2*iSeq+1], iL);
                    };
                };

                if (iSeq<nInd && ind1 == (int64)(ii-mapGen.genomeSAindexStart[iL])) {
                    SAi.writePacked(ii, indArray[2*iSeq]+iSeq);
                    for (uint ii0=ind0+1; ii0<ii; ii0++) {
                        SAi.writePacked(ii0, (indArray[2*iSeq]+iSeq) | mapGen.SAiMarkAbsentMaskC);
                    };
                    ++iSeq;
                    ind0=ii;
                } else {
                    iSeq=iSeq1;
                };
            } else {
                while (iSeq<nInd && indArray[2*iSeq]<iSA2) {
                    ++iSeq;
                };

                while (iSeq<nInd && indArray[2*iSeq]==iSA2) {
                    if (funCalcSAi(seq1[0]+indArray[2*iSeq+1], iL) >= (int64)(ii-mapGen.genomeSAindexStart[iL])) {
                        break;
                    };
                    ++iSeq;
                };

                SAi.writePacked(ii, iSA1+iSeq);

                for (uint ii0=ind0+1; ii0<ii; ii0++) {
                    SAi.writePacked(ii0, (iSA2+iSeq) | mapGen.SAiMarkAbsentMaskC);
                };
                ind0=ii;
            };
        };
    };

    for (uint64 iSeq=0; iSeq<nInd; iSeq++) {
        int64 ind1=0;
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {
            uint g=(uint) seq1[0][indArray[2*iSeq+1]+iL];
            ind1 <<= 2;
            if (g>3) {
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    ind1 += 3;
                    int64 ind2=mapGen.genomeSAindexStart[iL1]+ind1;
                    for (; ind2>=0; ind2--) {
                        if ((SAi[ind2] & mapGen.SAiMarkAbsentMaskC)==0) {
                            break;
                        };
                    };
                    SAi.writePacked(ind2, SAi[ind2] | mapGen.SAiMarkNmaskC);
                    ind1 <<= 2;
                };
                break;
            } else {
                ind1 += g;
            };
        };
    };

    time(&rawtime);
    P.inOut->logMain << timeMonthDayTime(rawtime) << "   Finished incremental SAi" << endl;
}
}

uint insertSeqSA(PackedArray & SA, PackedArray & SA1, PackedArray & SAi, char * G1, uint64 nG, uint64 nG1, uint64 nG2, Parameters & P, Genome &mapGen)
{//insert new sequences into the SA

    uint GstrandBit1 = (uint) floor(log(nG+nG1)/log(2))+1;
    if (GstrandBit1<32) GstrandBit1=32; //TODO: use simple access function for SA
    if ( GstrandBit1+1 > SA.wordLength)
    {//sequence is too long - GstrandBit changed
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: cannot insert sequence on the fly because of strand GstrandBit problem\n";
        errOut << "SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    uint N2bit= 1LLU << (SA.wordLength-1);

    char** seq1=new char*[2];

    #define GENOME_endFillL 16
    char* seqq=new char [4*nG1+3*GENOME_endFillL];//ends shouldbe filled with 5 to mark boundaries

    seq1[0]=seqq+GENOME_endFillL;//TODO: avoid defining an extra array, use reverse search
    seq1[1]=seqq+2*GENOME_endFillL+2*nG1;

    memset(seqq,GENOME_spacingChar,GENOME_endFillL);
    memset(seqq+2*nG1+GENOME_endFillL,GENOME_spacingChar,GENOME_endFillL);
    memset(seqq+4*nG1+2*GENOME_endFillL,GENOME_spacingChar,GENOME_endFillL);

    memcpy(seq1[0], G1, nG1);
    for (uint ii=0; ii<nG1; ii++)
    {//reverse complement sequence
        seq1[0][2*nG1-1-ii]=seq1[0][ii]<4 ? 3-seq1[0][ii] : seq1[0][ii];
    };
    complementSeqNumbers(seq1[0], seq1[1], 2*nG1);//complement

    uint64* indArray=new uint64[nG1*2*2+2];// for each base, 1st number - insertion place in SA, 2nd number - index, *2 for reverse compl
    uint64 nInd=0;//true number of new indices

    time_t rawtime;
    if (P.pGe.gInsertOverlayDeltaFile!="-") {
        nInd=readGenomeInsertDelta(P.pGe.gInsertOverlayDeltaFile, indArray, 2*nG1, nG, nG1, nG2, G1, SA, P, mapGen);
        time ( &rawtime );
        P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Loaded genome insert delta, number of new SA indices = "<<nInd<<endl;
    } else {

        #pragma omp parallel num_threads(P.runThreadN)
        #pragma omp for schedule (dynamic,64)
        for (uint ii=0; ii<2*nG1; ii++) {//find insertion points for each of the sequences

            if (seq1[0][ii]>3)
            {//no index for suffices starting with N
                indArray[ii*2]=-1;
            } else
            {
                indArray[ii*2] =  suffixArraySearch1(mapGen, seq1, ii, 10000, nG, (ii<nG1 ? true:false), 0, SA.length-1, 0, nG, nG1, nG2) ;
                indArray[ii*2+1] = ii;
            };
        };

        for (uint ii=0; ii<2*nG1; ii++) {//remove entries that cannot be inserted, this cannot be done in the parallel cycle above
            if (indArray[ii*2]!= (uint) -1) {
                indArray[nInd*2]=indArray[ii*2];
                indArray[nInd*2+1]=indArray[ii*2+1];
                ++nInd;
            };
        };

        time ( &rawtime );
        P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SA search, number of new SA indices = "<<nInd<<endl;

    /*//old-debug
    uint64* indArray1=new uint64[nG1*2*2+2];
    memcpy((void*) indArray1, (void*) indArray, 8*(nG1*2*2+2));
    g_funCompareUintAndSuffixes_G=seq1[0];
    qsort((void*) indArray1, nInd, 2*sizeof(uint64), funCompareUintAndSuffixes);
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished qsort - old " <<endl;
    */

        g_funCompareUintAndSuffixesMemcmp_G=seq1[0];
        g_funCompareUintAndSuffixesMemcmp_L=mapGen.pGe.gSuffixLengthMax/sizeof(uint64_t);
        qsort((void*) indArray, nInd, 2*sizeof(uint64_t), funCompareUintAndSuffixesMemcmp);

//     qsort((void*) indArray, nInd, 2*sizeof(uint64), funCompareUint2);
        time ( &rawtime );
        P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished qsort" <<endl;
        if (P.pGe.gInsertOutMode=="Delta") {
            const string deltaFileName=genomeInsertDeltaFilePath(P.pGe.gInsertOutDir);
            writeGenomeInsertDelta(deltaFileName, indArray, nInd, nG, nG1, nG2, G1, SA, P, mapGen);
            time ( &rawtime );
            P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Wrote genome insert delta: " << deltaFileName << endl;
            return nInd;
        };
    };



    /*//new sorting, 2-step: qsort for indArray, bucket sort for suffixes
    qsort((void*) indArray, nInd, 2*sizeof(uint64), funCompareUint2);
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished qsort"<<nInd<<endl;

    sortSuffixesBucket(seq1[0], (void*) indArray, nInd, 2*sizeof(uint64));
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished ordering suffixes"<<nInd<<endl;
    */

    /* //debug
    for (int ii=0;ii<2*nInd;ii++)
    {
        if (indArray[ii]!=indArray1[ii])
        {
            cout << ii <<"   "<< indArray[ii]  <<"   "<< indArray1[ii] <<endl;
        };
    };
    */

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished sorting SA indices"<<endl;

    if (P.pGe.gInsertOverlayDeltaFile!="-" && !P.pGe.gInsertOverlayGTFhasJunctions && P.runMode=="alignReads") {
        insertSeqSAi(SAi, seq1, indArray, nInd, P, mapGen);
        mapGen.genomeInsertSAsetup(indArray, nInd, nG, nG1, nG2, N2bit);
        mapGen.nSA=SA.length+nInd;

        time ( &rawtime );
        P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Using virtual genome insert SA overlay, expanded SA indices = " << mapGen.nSA << endl;
        return nInd;
    };

    indArray[2*nInd]=-999; //mark the last junction
    indArray[2*nInd+1]=-999; //mark the last junction

    SA1.defineBits(SA.wordLength,SA.length+nInd);

    /*testing
    PackedArray SAo;
    SAo.defineBits(mapGen.GstrandBit+1,mapGen.nSA+nInd);
    SAo.allocateArray();
    ifstream oldSAin("./DirTrue/SA");
    oldSAin.read(SAo.charArray,SAo.lengthByte);
    oldSAin.close();
    */

    vector<uint64> waveBoundaries;
    waveBoundaries.push_back(0);
    if (P.runThreadN>1 && SA.length>0 && SA.wordLength+7<=sizeof(uint)*8) {
        // SA points into the tail of SA1's larger buffer. Process waves from
        // left to right and buffer each wave before writing so expanded output
        // never overwrites unread source SA records.
        const uint64 waveTargetBytes=1LLU << 30;
        uint nWaves=(SA.length*sizeof(uint)+waveTargetBytes-1)/waveTargetBytes;
        if (nWaves<2) {
            nWaves=2;
        };
        if (nWaves>(uint) P.runThreadN) {
            nWaves=(uint) P.runThreadN;
        };
        if (nWaves>64) {
            nWaves=64;
        };
        for (uint ii=1; ii<nWaves; ii++) {
            uint64 oldTarget=SA.length*ii/nWaves;
            uint64 oldBoundary=genomeInsertFindAlignedBoundary(oldTarget, SA.length, indArray, nInd, SA.wordLength);
            if (oldBoundary>waveBoundaries.back() && oldBoundary<SA.length) {
                waveBoundaries.push_back(oldBoundary);
            };
        };
    };

    uint64 oldSerialStart=0;
    if (waveBoundaries.size()>1) {
        for (uint64 iWave=0; iWave<waveBoundaries.size()-1; iWave++) {
            uint64 oldWaveBegin=waveBoundaries[iWave];
            uint64 oldWaveEnd=waveBoundaries[iWave+1];
            vector<uint> oldSAvalues(oldWaveEnd-oldWaveBegin);

            #pragma omp parallel for schedule(static) num_threads(P.runThreadN)
            for (uint64 old1=oldWaveBegin; old1<oldWaveEnd; old1++) {
                oldSAvalues[old1-oldWaveBegin]=SA[old1];
            };

            vector<uint64> oldBoundaries;
            oldBoundaries.push_back(oldWaveBegin);
            for (uint ii=1; ii<(uint) P.runThreadN; ii++) {
                uint64 oldTarget=oldWaveBegin+(oldWaveEnd-oldWaveBegin)*ii/P.runThreadN;
                uint64 oldBoundary=genomeInsertFindAlignedBoundary(oldTarget, oldWaveEnd, indArray, nInd, SA.wordLength);
                if (oldBoundary>oldBoundaries.back() && oldBoundary<oldWaveEnd) {
                    oldBoundaries.push_back(oldBoundary);
                };
            };
            oldBoundaries.push_back(oldWaveEnd);

            #pragma omp parallel for schedule(static) num_threads(P.runThreadN)
            for (uint64 ii=0; ii<oldBoundaries.size()-1; ii++) {
                insertSeqSAwriteChunk(SA1, oldSAvalues, oldWaveBegin, SA.wordLength, indArray, nInd, oldBoundaries[ii], oldBoundaries[ii+1], nG, nG1, nG2, N2bit);
            };
        };
        oldSerialStart=waveBoundaries.back();
    };

    uint64 isa1=genomeInsertCountBefore(indArray, nInd, oldSerialStart);
    uint64 isa2=oldSerialStart+isa1;
    for (uint64 isa=oldSerialStart;isa<SA.length;isa++) {
        while (isa==indArray[isa1*2]) {//insert new index before the existing index
            SA1.writePacked(isa2,genomeInsertInsertedSAvalue(indArray[isa1*2+1], nG, nG1, nG2, N2bit));
            /*testing
            if (SA1[isa2]!=SAo[isa2]) {
               cout <<isa2 <<" "<< SA1[isa2]<<" "<<SAo[isa2]<<endl;
               //sleep(100);
            };
            */
            ++isa2; ++isa1;

        };

        SA1.writePacked(isa2,genomeInsertSAshift(SA[isa], nG, nG1, nG2, SA.wordLength-1)); //TODO make sure that the first sj index is not before the first array index
            /*testing
            if (SA1[isa2]!=SAo[isa2]) {
               cout <<isa2 <<" "<< SA1[isa2]<<" "<<SAo[isa2]<<endl;
               //sleep(100);
            };
            */
        ++isa2;
    };
    for (;isa1<nInd;isa1++)
    {//insert the last indices
        SA1.writePacked(isa2,genomeInsertInsertedSAvalue(indArray[isa1*2+1], nG, nG1, nG2, N2bit));
        ++isa2;
    };

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished inserting SA indices" << (waveBoundaries.size()>1 ? " in parallel" : "") <<endl;

    insertSeqSAi(SAi, seq1, indArray, nInd, P, mapGen);

    //change parameters, most parameters are already re-defined in sjdbPrepare.cpp
    SA.defineBits(mapGen.GstrandBit+1,SA.length+nInd);//same as SA2
    SA.pointArray(SA1.charArray);
    mapGen.nSA=SA.length;
    mapGen.nSAbyte=SA.lengthByte;

//     mapGen.sjGstart=mapGen.chrStart[mapGen.nChrReal];
//     memcpy(G+mapGen.chrStart[mapGen.nChrReal],seq1[0], nseq1[0]);


    return nInd;
};
