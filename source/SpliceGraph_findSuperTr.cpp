#include "SpliceGraph.h"
#include "sjAlignSplit.h"
#include "ReadAlign.h"

void SpliceGraph::findSuperTr(const char *readSeq, const char *readSeqRevCompl, const uint32 readLen, const string &readName, Genome &mapGen)
{//find the candidate superTranscripts: seed-and-rank algorithm implemented
		
	//hardcoded for now
	float seedCoverageThreshold = 0.02;
	float seedCoverageMinToMax = 0.5;
	uint32 seedMultMax = P.seedMultimapNmax;
    uint32 seedSpacing = 1;
	uint32 seedLen=mapGen.pGe.gSAindexNbases; //TODO: make user-definable	vector<uint32> seedSuperTr;
    
    vector<uint32> seedSuperTr;
	seedSuperTr.reserve(seedMultMax);
    memset(superTrSeedCount,0,sizeof(superTrSeedCount[0])*2*superTrome.N);
    
    for (uint32 iseed=0; iseed<readLen; iseed+=seedSpacing) {//loop through seeds
        //calculate index
        uint64 ind1=0;
        for (uint32 ii=iseed; ii<iseed+seedLen; ii++) {
            uint b=(uint64) readSeq[ii];
            if (b>3) {//N
                continue; //no mapping for seeds with Ns
            } else {
                ind1 <<= 2;
                ind1 += b;
            };
        };
        //find seed boundaries in SA
        uint64 iSA1=mapGen.SAi[mapGen.genomeSAindexStart[seedLen-1]+ind1]; // starting point for suffix array search.
        if ( (iSA1 & mapGen.SAiMarkAbsentMaskC) != 0) {//prefix does not exist, skip this seed
            continue;
        };
        uint64 iSA2;
        if (mapGen.genomeSAindexStart[seedLen-1]+ind1+1 < mapGen.genomeSAindexStart[seedLen]) {//we are not at the end of the SA
            iSA2=((mapGen.SAi[mapGen.genomeSAindexStart[seedLen-1]+ind1+1] & mapGen.SAiMarkNmask) & mapGen.SAiMarkAbsentMask) - 1;
        } else {
            iSA2=mapGen.nSA-1;
        };

		if (iSA2-iSA1>=seedMultMax) //this seed map too many times
			continue;
		
		seedSuperTr.clear();
		seedSuperTr.resize(iSA2-iSA1+1);
        for (uint64 isa=iSA1;isa<=iSA2;isa++) {//loop through seed SA boundaries
            uint64 a1 = mapGen.SA[isa];
            uint64 aStr = a1 >> mapGen.GstrandBit;
            a1 = a1 & mapGen.GstrandMask; //remove strand bit
            if (aStr==1)
                a1 = mapGen.nGenome - (seedLen+a1);
            
            if (a1>=mapGen.sjGstart) {//this is sj align. Probably do not need it...
                uint64 a1D, aLengthD, a1A, aLengthA, sj1;
                if ( sjAlignSplit(a1, seedLen, mapGen, a1D, aLengthD, a1A, aLengthA, sj1) ) {//if align crosses the junction, replace a1 with a1D
                    a1=a1D;
                } else {
                    continue;//this align does not cross the junction, no need to contine
                };
            };
            
            seedSuperTr[isa-iSA1]=(uint32)(aStr*superTrome.N  + mapGen.chrBin[a1 >> mapGen.pGe.gChrBinNbits]);
        };//loop through seed SA boundaries
		
		sort(seedSuperTr.begin(),seedSuperTr.end());
		uint32 su1prev=(uint32)-1;
		for (auto &su1 : seedSuperTr) {//this will eliminate multiple matches of the same seed into the same suTr
			if (su1!=su1prev) {
				superTrSeedCount[su1]++;
				su1prev=su1;
			};
		};
			
    };//loop through seeds
	
    //find max coverage
    typeSuperTrSeedCount countMax=0;
	//float countOverSuperTrLenMax=0;
    for (uint32 ii=0; ii<2*superTrome.N; ii++) {
        countMax=max(superTrSeedCount[ii], countMax);
		//countOverSuperTrLenMax=max(superTrSeedCount[ii]/float(superTr.length[ii%superTr.N]), countOverSuperTrLenMax);
    };
    
    if (countMax<readLen*seedCoverageThreshold/seedSpacing)//no good candidates, hard-coded for now
        return;
    
    uint32 nSuperTr=0;
    uint32 maxMaxScore=0;
    for (uint32 ii=0; ii<superTrome.N; ii++) {//selection cycle
				
        uint32 sutr1=ii%superTrome.N;
        uint32 str1=ii/superTrome.N;
		
        //if (superTrSeedCount[ii]<countOverSuperTrLenMax*superTr.length[sutr1]*seedCoverageMinToMax)
        //if (superTrSeedCount[ii]<countMax*seedCoverageMinToMax)
		if (superTrSeedCount[ii] < readLen*seedCoverageThreshold/seedSpacing || superTrSeedCount[ii]<countMax*seedCoverageMinToMax) {
			continue;
        };

		if (readLen>=100000 || superTrome.superTrs[sutr1].length*readLen>=1000000000) {//temporary restriction: will lift after redoing the scoringMatrix
			continue;		
        };
                
        //convert into trAll
        RA->trAll[nSuperTr]=RA->trArrayPointer+nSuperTr;
        RA->nWinTr[nSuperTr]=1;
        Transcript &trA = *RA->trAll[nSuperTr][0]; //transcript to fill
        trA=*RA->trInit;
        trA.Chr=sutr1;
        trA.Str=str1;
        
		uint32 swScore = 0;
        swScore = swScoreSpliced((str1==0 ? readSeq : readSeqRevCompl), readLen, superTrome.superTrs[sutr1], trA.cigar);
        
        //swTraceBack(alignInfo.aEnd, alignInfo.aStart);
		//float(superTrSeedCount[ii])/countMax
		//float(superTrSeedCount[ii])/superTr.length[sutr1]/countOverSuperTrLenMax <<'\t'<<

		cout << readName <<'\t'<< sutr1 <<'\t'<< str1 <<'\t'<< superTrome.superTrs[sutr1].length <<'\t'<< readLen <<'\t'<< float(superTrSeedCount[ii])/readLen*seedSpacing <<'\t'<<  
				float(superTrSeedCount[ii])/countMax <<'\t'<<
                swScore <<'\t'<< float(swScore)/readLen <<'\t'<< alignInfo.aStart[0] <<'\t'<< alignInfo.aEnd[0] <<'\t'<< alignInfo.aStart[1] <<'\t'<< alignInfo.aEnd[1] << endl;
        
        trA.maxScore=swScore;
        trA.nMatch=swScore;//TODO fix this, calculate number of matched bases
        trA.nExons=0;
        trA.gStart=mapGen.chrStart[trA.Chr]+alignInfo.aStart[1];
        trA.nMM=alignInfo.nMM;
        trA.lIns=alignInfo.nI;
        trA.lDel=alignInfo.nD;
        trA.rLength=alignInfo.nMap;
//         {//calculate blocks from rowCol and rowSJ
//             int32 iEx=-1;//current exon
//             for (uint32 row=alignInfo.aStart[0]; row<=alignInfo.aEnd[0]; row++) {
// 
//                 if (rowCol[row]>=0) {//this row has no mapped base (i.e. no block) = bases deleted from query
//                     if (iEx==-1 || rowCol[row]!=rowCol[row-1]+1) {//start new block
//                         ++iEx;
//                         trA.exons[iEx][EX_R]=row;
//                         trA.exons[iEx][EX_G]=mapGen.chrStart[trA.Chr]+rowCol[row];
//                         trA.exons[iEx][EX_L]=1;
//                         trA.canonSJ[iEx]=-1;
//                         trA.sjAnnot[iEx]=0;
//                     } else {
//                         trA.exons[iEx][EX_L]++;
//                     };
//                 };
//          
//                 //this needs to be done even if rowCol[row]==-1
//                 if (rowSJ[row][0]>=0) {//junction: create two 0-length blocks to encapsulate junction
//                     ++iEx;
//                     trA.exons[iEx][EX_R]=row+1;
//                     trA.exons[iEx][EX_G]=mapGen.chrStart[trA.Chr]+rowSJ[row][0]+1;
//                     trA.exons[iEx][EX_L]=0;
//                     trA.canonSJ[iEx]=1;//TODO actual motif
//                     trA.sjAnnot[iEx]=1;
//                     //if (rowSJ[row]!=rowCol[row+1]) {//create 0-length block to mark junction, on the next row
//                         ++iEx;
//                         trA.exons[iEx][EX_R]=row+1;
//                         trA.exons[iEx][EX_G]=mapGen.chrStart[trA.Chr]+rowSJ[row][1];
//                         trA.exons[iEx][EX_L]=0;
//                         trA.canonSJ[iEx]=-1;
//                         trA.sjAnnot[iEx]=0;
//                     //};
//                 };
//             };
//             
//             trA.nExons=iEx+1;
// 
//             for (uint32 iex=0; iex<trA.nExons; iex++) {
//                 trA.exons[iex][EX_iFrag]=0;
//                 trA.exons[iex][EX_sjA]=0;
//             };            
//         };       

        if (swScore>maxMaxScore) {
            maxMaxScore=swScore;
            RA->trBest=&trA;
        };
        nSuperTr++;
    };
    RA->nW=nSuperTr;
    return;
};

