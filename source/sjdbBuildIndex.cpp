#include "sjdbBuildIndex.h"
// #include "sjdbLoadFromStream.h"
// #include "sjdbPrepare.h"
#include "ErrorWarning.h"
#include "SuffixArrayFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include "IncludeDefine.h"
#include "streamFuns.h"
#include "binarySearch2.h"
#include "ErrorWarning.h"
#include <cmath>

#include "funCompareUintAndSuffixes.h"

namespace {

struct SAiPart {//entries [ii0,ii1) of SAi level iL; every part but a level's first starts at a present entry
    uint iL, ii0, ii1;
    bool first, last;//first or last part of its level
    uint V, k, lb;//not first: the first entry's value without marks, its offset in the level, the first junction with indArray[2*j]>=V
    uint iSJguess, iSJdry, ind0dry, orDry;//the dry run's starting state, final state and last present entry, OR of the values it would write
    uint iSJin, S, fill;//true state entering the part and after its first entry; the value the next part's first entry fills this part's trailing absent entries with
    vector<pair<uint,uint>> deferred;//writes whose 8 bytes reach a byte not wholly owned by this part
};

//One part of the serial SAi pass, entries [ii0,ii1) of level iL from state iSJ, ind0: the loop body is the serial one.
//write==false is a dry run: nothing is written, and orVals collects every value the serial loop would write.
//write==true: a write whose 8-byte read-modify-write touches a byte outside [byteLo,byteHi) goes to deferred;
//with fillTail, the absent entries after the last present one are filled with fillVal at the end.
uint saiScan(PackedArray &SAi, Genome &mapGen, const uint64 *indArray, char *Gsj, uint nInd, uint iL, uint ii0, uint ii1,
             uint iSJ, uint &ind0, bool write, uint byteLo, uint byteHi, bool fillTail, uint fillVal,
             vector<pair<uint,uint>> &deferred, uint &orVals)
{
    const uint start=mapGen.genomeSAindexStart[iL], wl=SAi.wordLength;
    auto put = [&] (uint ii, uint x) {
        const uint B=(ii*wl)>>3;
        if (B>=byteLo && B+8<=byteHi) {
            SAi.writePacked(ii,x);
        } else {
            deferred.push_back(make_pair(ii,x));
        };
    };
    for (uint ii=ii0; ii<ii1; ii++) {
        uint iSA1=SAi[ii];
        uint iSA2=iSA1 & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask;

        if ( iSJ<nInd && (iSA1 &  mapGen.SAiMarkAbsentMaskC)>0 )
        {//index missing from the old genome
            uint iSJ1=iSJ;
            int64 ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
            while (ind1 < (int64)(ii-start) && indArray[2*iSJ]-1<iSA2) {
                ++iSJ;
                ind1=funCalcSAi(Gsj+indArray[2*iSJ+1],iL);
            };
            if (ind1 == (int64)(ii-start) ) {
                const uint x=indArray[2*iSJ]-1+iSJ+1;
                if (write) {
                    put(ii,x);
                    for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                        put(ii0,x | mapGen.SAiMarkAbsentMaskC);
                    };
                } else {
                    orVals |= x | (ind0+1<ii ? mapGen.SAiMarkAbsentMaskC : 0);
                };
                ++iSJ;
                ind0=ii;
            } else {
                iSJ=iSJ1;
            };
        } else
        {//index was present in the old genome
            while (iSJ<nInd && indArray[2*iSJ]-1+1<iSA2) {//for this index insert "smaller" junctions
                ++iSJ;
            };

            while (iSJ<nInd && indArray[2*iSJ]-1+1==iSA2) {//special case, the index falls right behind SAi
                if (funCalcSAi(Gsj+indArray[2*iSJ+1],iL) >= (int64) (ii-start) ) {//this belongs to the next index
                    break;
                };
                ++iSJ;
            };

            if (write) {
                put(ii,iSA1+iSJ);
                for (uint ii0=ind0+1; ii0<ii; ii0++) {//fill all the absent indices with this value
                    put(ii0,(iSA2+iSJ) | mapGen.SAiMarkAbsentMaskC);
                };
            } else {
                orVals |= (iSA1+iSJ) | (ind0+1<ii ? (iSA2+iSJ) | mapGen.SAiMarkAbsentMaskC : 0);
            };
            ind0=ii;
        };
    };
    if (write && fillTail) {
        for (uint ii0=ind0+1; ii0<ii1; ii0++)
            put(ii0,fillVal);
    };
    return iSJ;
};

//The SAi pass of the junction insertion. Levels are independent; a level is split into parts that start at present entries.
//At a present entry of value V the serial loop moves iSJ past the junctions below V and the tied ones that belong to an
//earlier entry; from any state up to the first junction it stops at, found by a binary search, it reaches that junction.
//Each part's dry run starts from that guess; a serial check then gives every part its true starting state, re-running
//a part whose guess leads elsewhere, and the parts run again, writing.
//A value wider than an entry spills into the next entry, which the serial loop reads after the write; any such value
//sends the whole pass to the serial loop.
void sjdbInsertSAi(PackedArray &SAi, Genome &mapGen, const uint64 *indArray, char *Gsj, uint nInd, int nThreads, ostream &logMain)
{
    const uint nL=mapGen.pGe.gSAindexNbases, wl=SAi.wordLength;
    const uint saiPartMin=1LLU<<12;//entries

    vector<SAiPart> parts;
    bool split=false;
    if (nThreads>1) {
        const uint partN=max(saiPartMin, mapGen.genomeSAindexStart[nL]/(4*(uint)nThreads));
        for (uint iL=0; iL<nL; iL++) {
            const uint start=mapGen.genomeSAindexStart[iL], end=mapGen.genomeSAindexStart[iL+1];
            const uint nP=max((uint)1, (end-start)/partN);
            SAiPart p;
            p.iL=iL; p.ii0=start; p.first=true;
            parts.push_back(p);
            for (uint ip=1; ip<nP; ip++) {
                uint ii=max(start+ip*((end-start)/nP), parts.back().ii0+1);
                while (ii<end && (SAi[ii] & mapGen.SAiMarkAbsentMaskC)>0)
                    ++ii;
                if (ii>=end)
                    break;
                parts.back().ii1=ii;
                p.ii0=ii; p.first=false;
                parts.push_back(p);
                split=true;
            };
            parts.back().ii1=end;
        };
    };

    auto serialPass = [&] () {
        vector<pair<uint,uint>> none;
        uint orV=0;
        for (uint iL=0; iL<nL; iL++) {
            uint ind0=mapGen.genomeSAindexStart[iL]-1;//last index that was present in the old genome
            saiScan(SAi, mapGen, indArray, Gsj, nInd, iL, mapGen.genomeSAindexStart[iL], mapGen.genomeSAindexStart[iL+1],
                    0, ind0, true, 0, (uint) -1, false, 0, none, orV);
        };
    };

    if (!split) {
        serialPass();
        return;
    };

    //the serial loop's two while loops at a present entry, from state iSJ; indArray is sorted on its even elements
    auto presentLoops = [&] (const SAiPart &p, uint iSJ) -> uint {
        iSJ=max(iSJ,p.lb);
        while (iSJ<nInd && indArray[2*iSJ]==p.V) {
            if (funCalcSAi(Gsj+indArray[2*iSJ+1],p.iL) >= (int64) p.k)
                break;
            ++iSJ;
        };
        return iSJ;
    };

    for (auto &p : parts) {
        p.last = p.ii1==mapGen.genomeSAindexStart[p.iL+1];
        p.iSJguess=0;
        if (!p.first) {
            p.V=SAi[p.ii0] & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask;
            p.k=p.ii0-mapGen.genomeSAindexStart[p.iL];
            uint i1=0, i2=nInd;
            while (i1<i2) {
                const uint im=(i1+i2)/2;
                if (indArray[2*im]<p.V) {
                    i1=im+1;
                } else {
                    i2=im;
                };
            };
            p.lb=i1;
            p.iSJguess=presentLoops(p,p.lb);
        };
    };

    vector<uint> order(parts.size());//largest parts first
    for (uint ip=0; ip<parts.size(); ip++)
        order[ip]=ip;
    sort(order.begin(), order.end(), [&] (uint a, uint b) {return parts[a].ii1-parts[a].ii0 > parts[b].ii1-parts[b].ii0;});

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic,1)
    for (uint io=0; io<order.size(); io++) {
        SAiPart &p=parts[order[io]];
        p.ind0dry=p.ii0-1;
        p.orDry=0;
        p.iSJdry=saiScan(SAi, mapGen, indArray, Gsj, nInd, p.iL, p.ii0, p.ii1, p.iSJguess, p.ind0dry, false, 0, 0, false, 0, p.deferred, p.orDry);
    };

    uint orAll=0, nRerun=0;
    for (uint ip=0; ip<parts.size(); ip++) {
        SAiPart &p=parts[ip];
        if (p.first) {
            p.iSJin=0;
        } else {
            SAiPart &q=parts[ip-1];
            p.iSJin=q.iSJdry;
            p.S=presentLoops(p,p.iSJin);
            if (p.S!=p.iSJguess) {//the serial loop leaves another state at the first entry
                p.ind0dry=p.ii0-1;
                p.orDry=0;
                p.iSJdry=saiScan(SAi, mapGen, indArray, Gsj, nInd, p.iL, p.ii0, p.ii1, p.iSJin, p.ind0dry, false, 0, 0, false, 0, p.deferred, p.orDry);
                ++nRerun;
            };
            q.fill=(p.V+p.S) | mapGen.SAiMarkAbsentMaskC;
            if (q.ind0dry+1<q.ii1)
                orAll |= q.fill;
        };
        orAll |= p.orDry;
    };

    if ( (orAll>>wl)>0 ) {
        logMain << "   SAi: a value is wider than an entry, serial pass" <<endl;
        serialPass();
        return;
    };

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic,1)
    for (uint io=0; io<order.size(); io++) {
        SAiPart &p=parts[order[io]];
        uint ind0=p.ii0-1, orV=0;
        saiScan(SAi, mapGen, indArray, Gsj, nInd, p.iL, p.ii0, p.ii1, p.iSJin, ind0, true, (p.ii0*wl+7)>>3, (p.ii1*wl)>>3,
                !p.last, p.fill, p.deferred, orV);
    };
    for (auto &p : parts) {
        for (auto &d : p.deferred)
            SAi.writePacked(d.first,d.second);
    };
    logMain << "   SAi: " << parts.size() << " parts, " << nRerun << " re-run from their true starting state" <<endl;
};

//suffixArraySearch1(mapGen, s, S, N, -1LLU, true, 0, mapGen.nSA-1, 0) with the same midpoints, compares and result,
//except that a midpoint outside [iLo,iHi], the SA range of the suffix's longest present SAi prefix, is not compared:
//every entry before iLo is smaller than the suffix and every entry after iHi larger. A skipped midpoint sets its side's
//matched length to 0, so later compares may re-check bases the full search skips; their compRes is the same.
uint suffixArraySearchSAi(Genome &mapGen, PackedArray &SAi, char** s, uint S, uint N)
{
    const uint nSA=mapGen.nSA;
    uint iLo=0, iHi=nSA-1;

    const char *s0=s[0]+S;
    uint Lind=0, ind1=0;
    while (Lind<mapGen.pGe.gSAindexNbases && (unsigned char) s0[Lind]<4) {//prefix of A,C,G,T only
        ind1=(ind1<<2) + (unsigned char) s0[Lind];
        ++Lind;
    };
    while (Lind>0) {
        const uint iSAi=mapGen.genomeSAindexStart[Lind-1]+ind1;
        const uint v=SAi[iSAi];
        if ((v & mapGen.SAiMarkAbsentMaskC)==0) {//prefix present: its SA range starts at v and, when the next prefix is present, ends before that prefix's start, as in maxMappableLength2strands
            uint lo=v & mapGen.SAiMarkNmask, hi=nSA-1;
            if (iSAi+1 < mapGen.genomeSAindexStart[Lind]) {
                const uint w=SAi[iSAi+1];
                if ((w & mapGen.SAiMarkAbsentMaskC)==0) {
                    hi=(w & mapGen.SAiMarkNmask) - 1;
                } else if ((SAi[mapGen.genomeSAindexStart[Lind]-1] & mapGen.SAiMarkAbsentMaskC)==0) {//an absent next prefix holds the next present prefix's start, trusted only when the level ends present, since absent entries after a level's last present one can be stale; maxMappableLength2strands uses nSA-1 here
                    hi=(w & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask) - 1;
                };
            };
            if (lo<=hi && hi<nSA) {//otherwise the index is inconsistent: search the whole SA
                iLo=lo;
                iHi=hi;
            };
            break;
        };
        --Lind;
        ind1 >>= 2;
    };

    int compRes;
    uint i1=0, i2=nSA-1, L1=0, L2=0;
    if (i1>=iLo) {
        L1=compareSeqToGenome1(mapGen,s,S,N,0,i1,true,-1LLU,compRes);
        if (compRes<0)
            return 0;
    };
    if (i2<=iHi) {
        L2=compareSeqToGenome1(mapGen,s,S,N,0,i2,true,-1LLU,compRes);
        if (compRes>0)
            return -2llu;
    };
    uint L=min(L1,L2);
    while (i1+1<i2) {
        const uint i3=i1/2 + i2/2 + (i1%2 + i2%2)/2;//medianUint2
        if (i3<iLo) {
            i1=i3; L1=0;
        } else if (i3>iHi) {
            i2=i3; L2=0;
        } else {
            const uint L3=compareSeqToGenome1(mapGen,s,S,N,L,i3,true,-1LLU,compRes);
            if (L3==N)
                return i3;
            if (compRes>0) {
                i1=i3; L1=L3;
            } else if (compRes<0) {
                i2=i3; L2=L3;
            };
        };
        L=min(L1,L2);
    };
    return i2;
};

}//namespace

void sjdbBuildIndex (Parameters &P, char *Gsj, char *G, PackedArray &SA, PackedArray &SA2, PackedArray &SAi, Genome &mapGen, Genome &mapGen1) {

    #define SPACER_CHAR GENOME_spacingChar

    if (mapGen.sjdbN==0)
    {//no junctions to insert
        return;
    };

    time_t rawtime;
    time ( &rawtime );
    P.inOut->logMain   << timeMonthDayTime(rawtime) << " ..... inserting junctions into the genome indices" <<endl;
    *P.inOut->logStdOut  << timeMonthDayTime(rawtime) << " ..... inserting junctions into the genome indices" <<endl;

    uint nGsj=mapGen.sjdbLength*mapGen.sjdbN;
    for (uint ii=1; ii<=mapGen.sjdbN; ii++)
    {
        Gsj[ii*mapGen.sjdbLength-1]=SPACER_CHAR; //to make sure this is > than any genome char
    };
    Gsj[nGsj*2]=SPACER_CHAR;//mark the end of the text

    const int nThreads = P.runThreadN>1 ? P.runThreadN : 1;

    //Junction i is old junction i for every i (binarySearch2 finds it there), and SA keeps its word length and size:
    //the search would mark every suffix as old (sjNew==0, nInd==0) and SA2 would get SA's bytes unchanged (saUnchanged below),
    //so the reverse complement, the search and the sort are skipped. Gsj's second half is read only by them.
    bool noNewSJ = mapGen.sjdbN==mapGen1.sjdbN && mapGen.chrStart[mapGen.nChrReal]+nGsj==mapGen1.nGenome
                   && SA.wordLength==(uint)mapGen.GstrandBit+1 && SA.length==mapGen.nSA
                   && memcmp(mapGen.sjdbStart, mapGen1.sjdbStart, mapGen.sjdbN*sizeof(mapGen.sjdbStart[0]))==0
                   && memcmp(mapGen.sjdbEnd, mapGen1.sjdbEnd, mapGen.sjdbN*sizeof(mapGen.sjdbEnd[0]))==0;
    if (noNewSJ) {
        int sameInd=1;
        #pragma omp parallel for num_threads(nThreads) schedule(static) reduction(&&:sameInd)
        for (uint isj1=0; isj1<mapGen.sjdbN; isj1++)
            sameInd = sameInd && binarySearch2(mapGen.sjdbStart[isj1],mapGen.sjdbEnd[isj1],mapGen1.sjdbStart,mapGen1.sjdbEnd,mapGen1.sjdbN)==(int)isj1;
        noNewSJ = sameInd!=0;
    };

//     uint nIndicesSJ1=mapGen.sjdbOverhang;
    uint   nIndicesSJ1=mapGen.sjdbLength;//keep all indices - this is pre-2.4.1 of generating the genome

    char* G1c=NULL;
    uint32* oldSJind=NULL;
    uint64* indArray=NULL;
    uint64 sjNew=0;
    uint nInd=0;//true number of new indices

    if (!noNewSJ) {
        for (uint ii=0; ii<nGsj; ii++) {//reverse complement junction sequences
            Gsj[nGsj*2-1-ii]=Gsj[ii]<4 ? 3-Gsj[ii] : Gsj[ii]; //reverse complement
        };

        G1c=new char[nGsj*2+1];
        complementSeqNumbers(Gsj, G1c, nGsj*2+1);

        oldSJind=new uint32[mapGen1.sjdbN];

        indArray=new uint64[2*mapGen.sjdbN*(nIndicesSJ1+1)*2];//8+4 bytes for SA index and index in the genome * nJunction * nIndices per junction * 2 for reverse compl
        #pragma omp parallel num_threads(P.runThreadN)
        #pragma omp for schedule (dynamic,1000) reduction(+:sjNew)
        for (uint isj=0; isj<2*mapGen.sjdbN; isj++) {//find insertion points for each of the sequences

            char** seq1=new char*[2];
            seq1[0]=Gsj+isj*mapGen.sjdbLength;
            seq1[1]=G1c+isj*mapGen.sjdbLength;

            uint isj1=isj<mapGen.sjdbN ? isj : 2*mapGen.sjdbN-1-isj;
            int sjdbInd = mapGen1.sjdbN==0 ? -1 : binarySearch2(mapGen.sjdbStart[isj1],mapGen.sjdbEnd[isj1],mapGen1.sjdbStart,mapGen1.sjdbEnd,mapGen1.sjdbN);
            if (sjdbInd<0)
            {//count new junctions
                ++sjNew;
            } else
            {//record new index of the old junctions
                oldSJind[sjdbInd]=isj1;
            };

            for (uint istart1=0; istart1<nIndicesSJ1;istart1++) {

                uint istart=istart1;
//             uint istart=isj<mapGen.sjdbN ? istart1 : istart1+1; //for rev-compl junction, shift by one base to start with the 1st non-spacer base
                uint ind1=2*(isj*nIndicesSJ1+istart1);
                if (sjdbInd>=0 || seq1[0][istart]>3)
                {//no index for already included junctions, or suffices starting with N
                    indArray[ind1]=-1;
                } else
                {
                    //indArray[ind1] =  suffixArraySearch(seq1, istart, mapGen.sjdbLength-istart1, G, SA, true, 0, mapGen.nSA-1, 0, P) ;
                    indArray[ind1] =  suffixArraySearchSAi(mapGen, SAi, seq1, istart, 10000) ;
                    //-1LLU results in suffixes for the new junctions to be always included in SA *after* the suffixes of the old junctions
                    //for identical suffixes, this may result in unstable ordering
                    indArray[ind1+1] = isj*mapGen.sjdbLength+istart;
                };
            };
        };
        sjNew = sjNew/2;//novel junctions were double counted on two strands
    };

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SA search: number of new junctions=" << sjNew <<", old junctions="<<mapGen.sjdbN-sjNew<<endl;

    g_funCompareUintAndSuffixes_G=Gsj;
    if (!noNewSJ) {
        for (uint ii=0; ii<2*mapGen.sjdbN*nIndicesSJ1; ii++) {//remove entries that cannot be inserted, this cannot be done in the parallel cycle above
            if (indArray[ii*2]!= (uint) -1) {
                indArray[nInd*2]=indArray[ii*2];
                indArray[nInd*2+1]=indArray[ii*2+1];
                ++nInd;
            };
        };

        qsort((void*) indArray, nInd, 2*sizeof(uint64), funCompareUintAndSuffixes);
    };
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished sorting SA indicesL nInd="<<nInd <<endl;

    if (!noNewSJ) {
        indArray[2*nInd]=-999; //mark the last junction
        indArray[2*nInd+1]=-999; //mark the last junction
    };

    mapGen.nGenome=mapGen.chrStart[mapGen.nChrReal]+nGsj;
    mapGen.nSA+=nInd;

    char GstrandBit1 = (char) (uint) floor(log(mapGen.nGenome)/log(2))+1;
    if (GstrandBit1<32) GstrandBit1=32; //TODO: use simple access function for SA

    P.inOut->logMain <<"Genome size with junctions="<<mapGen.nGenome<<"  "<<mapGen.chrStart[mapGen.nChrReal]<<"   "<<nGsj<<"\n";
    P.inOut->logMain <<"GstrandBit1="<<int(GstrandBit1)<<"   GstrandBit="<<int(mapGen.GstrandBit)<<"\n";

    if ( GstrandBit1 > mapGen.GstrandBit)
    {//too many junctions were added - GstrandBit changed
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem\n";
        errOut << "SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    SA2.defineBits(mapGen.GstrandBit+1,mapGen.nSA);
    uint nGsjNew=sjNew*mapGen.sjdbLength; //this is the actual number of bytes added to the genome, while nGsj is the total size of all junctions

    uint N2bit= 1LLU << mapGen.GstrandBit;
    uint strandMask=~N2bit;

    /*testing
    PackedArray SAo;
    SAo.defineBits(mapGen.GstrandBit+1,mapGen.nSA);
    SAo.allocateArray();
    ifstream oldSAin("./DirTrue/SA");
    oldSAin.read(SAo.charArray,SAo.lengthByte);
    oldSAin.close();
    */

    const uint nSA1=mapGen1.nSA, nGenome1=mapGen1.nGenome, nGenomeNew=mapGen.nGenome;
    const uint sjGstart=mapGen.chrStart[mapGen.nChrReal], sjdbLength=mapGen.sjdbLength;

    auto newIndexSJ = [&] (uint isj) -> uint {//SA value of the inserted junction suffix isj
        uint ind1=indArray[isj*2+1];
        if (ind1<nGsj) {
            ind1+=sjGstart;
        } else {//reverse strand
            ind1=(ind1-nGsj) | N2bit;
        };
        return ind1;
    };

    auto newIndexOld = [&] (uint ind1) -> uint {//new SA value of the old SA value ind1
        if ( (ind1 & N2bit)>0 )
        {//- strand
            uint ind1s = nGenome1 - (ind1 & strandMask);
            if (ind1s>=sjGstart)
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1s-sjGstart)/sjdbLength;//old junction index
                ind1s += (oldSJind[sj1]-sj1)*sjdbLength;
                ind1 = (nGenomeNew - ind1s) | N2bit;
            } else
            {
                ind1+=nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
        } else
        {//+ strand
            if (ind1>=sjGstart)
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1-sjGstart)/sjdbLength;//old junction index
                ind1 += (oldSJind[sj1]-sj1)*sjdbLength;
            };
        };
        return ind1;
    };

    auto nSJbefore = [&] (uint isa) -> uint {//number of junction suffixes inserted before old index isa
        uint i1=0, i2=nInd;
        while (i1<i2) {
            uint im=(i1+i2)/2;
            if (indArray[2*im]<isa) {
                i1=im+1;
            } else {
                i2=im;
            };
        };
        return i1;
    };

    //No suffix inserted and every old junction kept its index (re-inserting the index's own junctions):
    //newIndexOld is the identity, so SA2 holds SA's bytes unchanged.
    bool saUnchanged = noNewSJ;
    if (!noNewSJ) {
        saUnchanged = nInd==0 && nGsjNew==0 && nGenomeNew==nGenome1 && mapGen.sjdbN==mapGen1.sjdbN
                      && SA.wordLength==SA2.wordLength && SA.length==SA2.length;
        for (uint isj1=0; saUnchanged && isj1<mapGen.sjdbN; isj1++)
            saUnchanged = binarySearch2(mapGen.sjdbStart[isj1],mapGen.sjdbEnd[isj1],mapGen1.sjdbStart,mapGen1.sjdbEnd,mapGen1.sjdbN)==(int)isj1;
    };

    if (saUnchanged) {
        //SA2 starts gapB bytes below SA in the same buffer: a round copies at most gapB bytes, so it writes only bytes already copied
        char *saFrom=SA.charArray, *saTo=SA2.charArray;
        const uint nB=SA.lengthByte, gapB = saFrom>saTo ? (uint) (saFrom-saTo) : 0;
        if (nThreads==1 || gapB<(1LLU<<26)) {
            memmove(saTo,saFrom,nB);
        } else {
            for (uint r0=0; r0<nB; r0+=gapB) {
                const uint r1=min(nB,r0+gapB), partB=(r1-r0+nThreads-1)/nThreads;
                #pragma omp parallel for num_threads(nThreads) schedule(static,1)
                for (int it=0; it<nThreads; it++) {
                    const uint c0=min(r1,r0+it*partB), c1=min(r1,c0+partB);
                    if (c1>c0)
                        memcpy(saTo+c0,saFrom+c0,c1-c0);
                };
            };
        };
        P.inOut->logMain << "   No suffixes inserted, old junctions keep their indices: SA moved unchanged" <<endl;
    } else {
        //SA2 starts below SA in the same buffer, so writing SA2 overwrites SA entries that were already read.
        //Rounds: each thread reads its block of SA, barrier, each thread writes its part of SA2, barrier.
        //A part's entries are packed in a register and stored as whole 64-bit words of SA2 (little-endian, as PackedArray);
        //its first and last words may hold bits of the neighbouring parts, so one thread merges those after the barrier.
        const uint blockN=1LLU<<18;//old SA entries per thread per round
        vector<uint> edgeWord(2*nThreads), edgeVal(2*nThreads), edgeMask(2*nThreads), edgeN(nThreads);

        #pragma omp parallel num_threads(nThreads)
        {
            const int nTeam=omp_get_num_threads();//OpenMP may start fewer threads than requested
            const int iThread=omp_get_thread_num();
            const uint roundN=blockN*nTeam;
            //private copies: the stores below may alias shared variables, which would then be re-read for every entry
            const char *saC=SA.charArray;
            const uint saWL=SA.wordLength, saMask=(~0LLU)>>(64-saWL);
            const uint gBit=mapGen.GstrandBit, sMask=strandMask, nG1=nGenome1, sjG=sjGstart, shiftMinus=nGsjNew;
            char *sa2C=SA2.charArray;
            const uint wl=SA2.wordLength;
            const uint64 *indA=indArray;
            vector<uint> newSA(blockN);
            uint *newSAp=newSA.data();
            for (uint isaRound=0; isaRound<nSA1; isaRound+=roundN) {
                const uint isaStart=min(nSA1, isaRound+iThread*blockN);
                const uint isaEnd=min(nSA1, isaStart+blockN);
                uint b=isaStart*saWL;
                for (uint isa=isaStart; isa<isaEnd; isa++, b+=saWL) {
                    uint ind1;
                    memcpy(&ind1, saC+(b>>3), sizeof(ind1));
                    ind1=(ind1>>(b&7)) & saMask;//SA[isa]
                    //newIndexOld(ind1), selecting the strand by a mask: the strand bit of consecutive entries is unpredictable.
                    //A suffix starting in the old junctions (rare) takes newIndexOld itself.
                    const uint minusMask=0-((ind1>>gBit)&1);
                    const uint gPos=(ind1 & ~minusMask) | ((nG1-(ind1 & sMask)) & minusMask);
                    newSAp[isa-isaStart] = gPos<sjG ? ind1+(shiftMinus & minusMask) : newIndexOld(ind1);
                };

                #pragma omp barrier
                uint isj=nSJbefore(isaStart);
                uint nextIns=indA[2*isj];
                //acc holds the bits of SA2 word w below bit accN that belong to this part; w0 is the part's first word
                const uint bit0=(isaStart+isj)*wl, w0=bit0>>6, lo0=bit0&63;
                uint w=w0, accN=lo0, acc=0, nEdge=0;
                auto put = [&] (uint ind1) {
                    acc |= ind1<<accN;
                    accN += wl;
                    if (accN>=64) {//word w is complete
                        if (w==w0) {
                            edgeWord[2*iThread]=w; edgeVal[2*iThread]=acc; edgeMask[2*iThread]=(~0LLU)<<lo0;
                            nEdge=1;
                        } else {
                            memcpy(sa2C+8*w, &acc, sizeof(acc));
                        };
                        ++w;
                        accN -= 64;
                        acc = ind1>>(wl-accN);//the bits that did not fit, 0 when accN==0
                    };
                };
                uint isa=isaStart;
                while (isa<isaEnd) {
                    while (isa==nextIns) {//insert sj index before the existing index
                        put(newIndexSJ(isj));
                        ++isj;
                        nextIns=indA[2*isj];
                    };
                    const uint runEnd=min(isaEnd, nextIns);
                    for (; isa<runEnd; isa++)
                        put(newSAp[isa-isaStart]);
                };
                if (w==w0) {//the part's bits lie in its first word, [lo0, accN)
                    if (accN>lo0) {
                        edgeWord[2*iThread]=w; edgeVal[2*iThread]=acc; edgeMask[2*iThread]=((~0LLU)>>(64-accN)) & ((~0LLU)<<lo0);
                        nEdge=1;
                    };
                } else if (accN>0) {//last word, partly filled
                    edgeWord[2*iThread+1]=w; edgeVal[2*iThread+1]=acc; edgeMask[2*iThread+1]=(~0LLU)>>(64-accN);
                    nEdge=2;
                };
                edgeN[iThread]=nEdge;

                #pragma omp barrier
                #pragma omp single
                {
                    for (int it=0; it<nTeam; it++) {
                        for (uint ie=0; ie<edgeN[it]; ie++) {
                            const uint ew=edgeWord[2*it+ie], ev=edgeVal[2*it+ie], em=edgeMask[2*it+ie];
                            for (uint k=0; k<8; k++) {//only the bytes holding this part's bits: a byte past the last entry may lie past the array
                                const unsigned char mb=(unsigned char) (em>>(8*k));
                                if (mb!=0) {
                                    unsigned char *c=(unsigned char*) sa2C+8*ew+k;
                                    *c=(unsigned char) ((*c & ~mb) | ((ev>>(8*k)) & mb));
                                };
                            };
                        };
                    };
                };
            };
        };

        uint isj=nSJbefore(nSA1), isa2=nSA1+isj;
        for (;isj<nInd;isj++) {//insert last new indices after the last old index
            SA2.writePacked(isa2,newIndexSJ(isj));
            ++isa2;
        };
    };

    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished inserting junction indices" <<endl;

    //SAi insertions; with nInd==0 every entry would be rewritten with its own value
    if (nInd>0)
        sjdbInsertSAi(SAi, mapGen, indArray, Gsj, nInd, nThreads, P.inOut->logMain);

    for (uint isj=0;isj<nInd;isj++) {
        int64 ind1=0;
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {
            uint g=(uint) Gsj[indArray[2*isj+1]+iL];
            ind1 <<= 2;
            if (g>3) {//this iSA contains N, need to mark the previous
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    ind1+=3;
                    int64 ind2=mapGen.genomeSAindexStart[iL1]+ind1;
                    for (; ind2>=0; ind2--) {//find previous index that is not absent
                        if ( (SAi[ind2] & mapGen.SAiMarkAbsentMaskC)==0 ) {
                            break;
                        };
                    };
                    SAi.writePacked(ind2,SAi[ind2] | mapGen.SAiMarkNmaskC);
                    ind1 <<= 2;
                };
                break;
            } else {
                ind1 += g;
            };
        };
    };
    time ( &rawtime );
    P.inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished SAi" <<endl;

    //change parameters, most parameters are already re-defined in sjdbPrepare.cpp
    SA.defineBits(mapGen.GstrandBit+1,mapGen.nSA);//same as SA2
    SA.pointArray(SA2.charArray);
    mapGen.nSAbyte=SA.lengthByte;
    mapGen.sjGstart=mapGen.chrStart[mapGen.nChrReal];
    memcpy(G+mapGen.chrStart[mapGen.nChrReal],Gsj, nGsj);

    /* testing
    PackedArray SAio=SAi;
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(P.pGe.gSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();


    for (uint ii=0;ii<mapGen.nSA;ii++) {
        if (SA2[ii]!=SAo[ii]) {
            cout <<ii <<" "<< SA2[ii]<<" "<<SAo[ii]<<endl;
        };
    };


    for (uint iL=0; iL < P.pGe.gSAindexNbases; iL++) {
        for (uint ii=mapGen.genomeSAindexStart[iL];ii<mapGen.genomeSAindexStart[iL+1]; ii++) {//scan through the longets index
                if ( SAio[ii]!=SAi[ii] ) {
                    cout <<ii<<" "<<SAio[ii]<<" "<<SAi[ii]<<endl;
                };
        };
    };
    */

    /*
    ofstream genomeOut("/home/dobin/Genome");
    fstreamWriteBig(genomeOut,G,mapGen.nGenome+nGsj,"777","777",P);
    genomeOut.close();
    genomeOut.open("/home/dobin/SA");
    fstreamWriteBig(genomeOut,SA2.charArray,SA2.lengthByte,"777","777",P);
    genomeOut.close();
    */

    delete [] indArray;
    delete [] G1c;
    delete [] oldSJind;

};
