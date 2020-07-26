#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include <unordered_map>
#include "SoloCommon.h"
#include <bitset>

#define def_MarkNoColor  (uint32) -1

inline int funCompareSolo1 (const void *a, const void *b) {
    uint32 *va= (uint32*) a;
    uint32 *vb= (uint32*) b;

    if (va[1]>vb[1]) {
		return 1;
	} else if (va[1]<vb[1]) {
		return -1;
	} else if (va[0]>vb[0]){
		return 1;
	} else if (va[0]<vb[0]){
		return -1;
	} else {
		return 0;
	};
};

void collapseUMIwith1MMlowHalf(uint32 *umiArr, uint32 umiArrayStride, uint32 umiMaskLow, uint32 nU0, uint32 &nU1, uint32 &nU2, uint32 &nC, vector<array<uint32,2>> &vC)
{
    const uint32 bitTop=1<<31;
    const uint32 bitTopMask=~bitTop;

    for (uint32 iu=0; iu<umiArrayStride*nU0; iu+=umiArrayStride) {//each UMI
        uint32 iuu=iu+umiArrayStride;
        for (; iuu<umiArrayStride*nU0; iuu+=umiArrayStride) {//compare to all UMIs down

            uint32 uuXor=umiArr[iu] ^ umiArr[iuu];

            if ( uuXor > umiMaskLow)
                break; //upper half is different

            if (uuXor >> (__builtin_ctz(uuXor)/2)*2 > 3) //shift by even number of trailing zeros
                continue;//>1MM

            //1MM UMI

            //graph coloring
            if ( umiArr[iu+2] == def_MarkNoColor && umiArr[iuu+2] == def_MarkNoColor ) {//no color
                //new color
                umiArr[iu+2] = nC;
                umiArr[iuu+2] = nC;
                ++nC;
                nU1 -= 2;//subtract the duplicated UMIs
            } else if ( umiArr[iu+2] == def_MarkNoColor ) {
                umiArr[iu+2] = umiArr[iuu+2];
                --nU1;//subtract the duplicated UMIs
            } else if ( umiArr[iuu+2] == def_MarkNoColor ) {
                umiArr[iuu+2] = umiArr[iu+2];
                --nU1;//subtract the duplicated UMIs
            } else {//both color
                if (umiArr[iuu+2] != umiArr[iu+2]) {//color conflict
                    //uint32 p[2]={umiArr[iu+2],umiArr[iuu+2]};
                    vC.push_back({umiArr[iu+2],umiArr[iuu+2]});
                    //vC.push_back({umiArr[iuu+2],umiArr[iu+2]});
                };
            };

            //directional collapse
            if ( (umiArr[iuu+1] & bitTop) == 0 && (umiArr[iu+1] & bitTopMask)>(2*(umiArr[iuu+1] & bitTopMask)-1) ) {//iuu is duplicate of iu
                umiArr[iuu+1] |= bitTop;
                --nU2;//subtract the duplicated UMIs
            } else if ( (umiArr[iu+1] & bitTop) == 0 && (umiArr[iuu+1] & bitTopMask)>(2*(umiArr[iu+1] & bitTopMask)-1) ) {//iu is duplicate of iuu
                umiArr[iu+1] |= bitTop;
                --nU2;//subtract the duplicated UMIs
            };
        };
    };
};

void graphDepthFirstSearch(uint32 n, vector<vector<uint32>> &nodeEdges, vector <uint32> &nodeColor) 
{
    for (const auto &nn : nodeEdges[n]) {
        if (nodeColor[nn]==(uint32)-1) {//node not visited
            nodeColor[nn]=nodeColor[n];
            graphDepthFirstSearch(nn,nodeEdges,nodeColor);
        };
    };
};

uint32 graphNumberOfConnectedComponents(uint32 N, vector<array<uint32,2>> V, vector<uint32> &nodeColor) 
{//find number of connected components
    //N=number of nodes
    //V=edges, list of connected nodes, each pair of nodes listed once
    //simple recursive DFS
   
    nodeColor.resize(N,(uint32)-1); //new color (connected component) for each node (each original color)

    if (V.size()==0)
        return N;

    vector<vector<uint32>> nodeEdges (N);
    for (uint32 ii=0; ii<V.size(); ii++) {
        nodeEdges[V[ii][0]].push_back(V[ii][1]);
        nodeEdges[V[ii][1]].push_back(V[ii][0]);
    };
    
    uint32 nConnComp=0;
    for (uint32 ii=0; ii<N; ii++) {
        if (nodeEdges[ii].size()==0) {//this node is not connected, no need to check. Save time because this happens often
            ++nConnComp;
        } else if (nodeColor[ii]==(uint32)-1) {//node not visited
            ++nConnComp;
            nodeColor[ii]=ii;
            graphDepthFirstSearch(ii,nodeEdges,nodeColor);
        };
    };
    return nConnComp;
};

void SoloFeature::collapseUMI(uint32 iCB, uint32 *umiArray) 
{
                                 
    uint32 *rGU=rCBp[iCB];
    uint32 rN=nReadPerCB[iCB];

    if (featureType==SoloFeatureTypes::Velocyto) {//special processing, by UMI
        qsort(rGU, rN, rguStride*sizeof(uint32), funCompareTypeShift<uint32, 1>);//sort by the UMI
        
//         unordered_map<uint32, array<uint32,3>>> geneCountHash;
        
        uint32 g1=(uint32)-1;//gene undefined for the first entry, for correct multigene identification
        bool multigene=false;//initial condition, to process first record
        uint32 gVelo=0;
        
        uint32 iRrec=0;
        for (uint32 iR=0; iR<rN*rguStride; iR+=rguStride) {//cycle over reads
            
            //process UMI
            uint32 g2=((rGU[iR] << velocytoTypeGeneBits) >> velocytoTypeGeneBits);//remove top 2 bits
            gVelo = gVelo | ( rGU[iR] >> velocytoTypeGeneBitShift );
            if (g2!=g1 && g1+1!=0)
                multigene=true;
            g1=g2;
            
            //IF end of array OR new UMI in the next record
            if (iR==(rN-1)*rguStride || rGU[iR+1]!=rGU[iR+3]) {//output this gene IF: past the end of array OR peeked ahead, and the umi is different
                if (!multigene) {//record this gene
                    rGU[iRrec]=g1;
                    ++iRrec;

                    bitset<velocytoTypeGeneBits> gV (gVelo);
                    if (!gV.test(AlignVsTranscript::ExonIntronSpan)) {//all UMIs are spanning models
                        rGU[iRrec]=2; //unspliced 
                    } else if (gV.test(AlignVsTranscript::Concordant)) {//>=1 purely exonic tr 
                        if (!gV.test(AlignVsTranscript::Intron) && !gV.test(AlignVsTranscript::ExonIntron)) {//0 purely intronic && 0 mixed
                            rGU[iRrec]=1; //spliced 
                        } else {//>=1 purely exonic and >=1 purely intronic or mixed
                            rGU[iRrec]=3; //ambiguous
                        };
                    } else {//0 exonic, >=1 intronic and/or >=1 mixed
                        rGU[iRrec]=2;//unspliced
                    };
                    ++iRrec;
                };

                g1=(uint32)-1;
                gVelo=0;
                multigene=false;
            };
        };
        
        countCellGeneUMIindex[iCB+1]=countCellGeneUMIindex[iCB];
        nGenePerCB[iCB]=0;
        nUMIperCB[iCB]=iRrec/2;
        
        if (iRrec==0) //all multigenes. nothing to record
            return;
    
        if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB] + (iRrec/2)*countMatStride) //allocated vector too small
            countCellGeneUMI.resize(countCellGeneUMI.size()*2);        
        
        qsort(rGU, iRrec/2, 2*sizeof(uint32), funCompareNumbers<uint32>);//sort by the gene

        countCellGeneUMIindex[iCB+1] -= countMatStride;//it will be added back for the first entry
        g1=(uint32)-1;//for the first entry to trigger new gene record
        for (uint32 iR=0; iR<iRrec; iR+=2) {
            if (rGU[iR]!=g1) {//new gene
                g1=rGU[iR];
                nGenePerCB[iCB]++;
                
                countCellGeneUMIindex[iCB+1] += countMatStride;//iCB+1 accumulates the index                
                countCellGeneUMI[countCellGeneUMIindex[iCB+1]+0]=g1;
            };
            countCellGeneUMI[countCellGeneUMIindex[iCB+1]+rGU[iR+1]]++;
        };
        countCellGeneUMIindex[iCB+1] += countMatStride;
        return;
    };    
    
    qsort(rGU,rN,rguStride*sizeof(uint32),funCompareNumbers<uint32>); //sort by gene index

    //compact reads per gene
    uint32 gid1=-1;//current gID
    uint32 nGenes=0; //number of genes
    uint32 *gID = new uint32[min(featuresNumber,rN)+1]; //gene IDs
    uint32 *gReadS = new uint32[min(featuresNumber,rN)+1]; //start of gene reads TODO: allocate this array in the 2nd half of rGU
    for (uint32 iR=0; iR<rN*rguStride; iR+=rguStride) {
        if (rGU[iR+rguG]!=gid1) {//record gene boundary
            gReadS[nGenes]=iR;
            gid1=rGU[iR+rguG];
            gID[nGenes]=gid1;
            ++nGenes;
        };
    };
    gReadS[nGenes]=rguStride*rN;//so that gReadS[nGenes]-gReadS[nGenes-1] is the number of reads for nGenes, see below in qsort

    unordered_map <uint32, unordered_map<uint32,uint32>> umiGeneHash, umiGeneHash1;
                   //UMI                 //Gene //Count
	unordered_map <uint32,uint32> geneCounts;

    /*
	if (pSolo.umiFiltering.MultiGeneUMI) {
		for (uint32 iR=0; iR<rN*rguStride; iR+=rguStride) {
			umiGeneHash[rGU[iR+1]][rGU[iR]]++;
		};
	};
    */

    if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB] + nGenes*countMatStride) //allocated vector too small
        countCellGeneUMI.resize(countCellGeneUMI.size()*2);
    
    nGenePerCB[iCB]=0;
    nUMIperCB[iCB]=0;
    countCellGeneUMIindex[iCB+1]=countCellGeneUMIindex[iCB];
    /////////////////////////////////////////////
    /////////// main cycle over genes
    for (uint32 iG=0; iG<nGenes; iG++) {//collapse UMIs for each gene
        uint32 *rGU1=rGU+gReadS[iG];

        if (gReadS[iG+1]==gReadS[iG])
        	continue; //no reads

        qsort(rGU1, (gReadS[iG+1]-gReadS[iG])/rguStride, rguStride*sizeof(uint32), funCompareTypeShift<uint32,rguU>);
        
        //exact collapse
        uint32 iR1=-umiArrayStride; //number of distinct UMIs for this gene
        uint32 u1=-1;
        for (uint32 iR=rguU; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//count and collapse identical UMIs
            
            if (rGU1[iR]!=u1) {
                iR1 += umiArrayStride;
                u1=rGU1[iR];
                umiArray[iR1]=u1;
                umiArray[iR1+1]=0;
                //umiArray[iR1+2]=def_MarkNoColor; //marks no color for graph
                //umiArray[iR1+2]=u1;
            };
            umiArray[iR1+1]++;
            //if ( umiArray[iR1+1]>nRumiMax) nRumiMax=umiArray[iR1+1];
        };

        uint32 nU0=(iR1+umiArrayStride)/umiArrayStride;
        uint32 nU1=nU0, nU2=nU0;//2 types of 1MM collapsing

        //testing exact matching CR counts.
		qsort(umiArray, nU0, umiArrayStride*sizeof(uint32), funCompareSolo1);
		for (uint64 iu=0; iu<(nU0-1)*umiArrayStride; iu+=umiArrayStride) {

			uint64 iuu;
			for (iuu=(nU0-1)*umiArrayStride; iuu>iu; iuu-=umiArrayStride) {

				uint32 uuXor=umiArray[iu+0] ^ umiArray[iuu+0];

				if ( (uuXor >> (__builtin_ctz(uuXor)/2)*2) <= 3 ) {//shift by even number of trailing zeros
					umiArray[iu+0]=umiArray[iuu+0]; //1MM
					break; //1MM
				};
			};
		};

	    if (pSolo.umiFiltering.MultiGeneUMI) {
	    	for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
	    		umiGeneHash[umiArray[iu+0]][gID[iG]]+=umiArray[iu+1];
	    	};
	    } else {
			qsort(umiArray, nU0, umiArrayStride*sizeof(uint32), funCompareNumbers<uint32>);
			nU1=1;
			for (uint64 iu=umiArrayStride; iu<nU0*umiArrayStride; iu+=umiArrayStride) {
				if (umiArray[iu+0]!=umiArray[iu+0-umiArrayStride]) {
					nU1++;
				};
			};
			nGenePerCB[iCB]++;
			nUMIperCB[iCB]+=nU1;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+0]=gID[iG];
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+1]=nU0;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+2]=nU1;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+3]=nU2;
			countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;//iCB+1 accumulates the index
	    };
    };

    if (pSolo.umiFiltering.MultiGeneUMI) {
		for (auto &iu: umiGeneHash) {
			uint32 maxu=0, maxg=-1;
			for (auto &ig : iu.second) {
				if (ig.second>maxu) {
					maxu=ig.second;
					maxg=ig.first;
				} else if (ig.second==maxu) {
					maxg=-1;
				};
			};
			if (maxg+1!=0)
				geneCounts[maxg]++;
		};

		for (auto &ig: geneCounts) {
			nGenePerCB[iCB]++;
			nUMIperCB[iCB]+=ig.second;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+0]=ig.first;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+1]=0;//probably won't need to compute multiple counts at once
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+2]=ig.second;
			countCellGeneUMI[countCellGeneUMIindex[iCB+1]+3]=0;
			countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB+1] + countMatStride;//iCB+1 accumulates the index
		};

    };
    /*if (readInfo.size()>0 && nU0<=1) {//fill in readInfo: CB,UMI
        for (uint32 iR=0; iR<gReadS[iG+1]-gReadS[iG]; iR+=rguStride) {//cycle over reads
            if (rGU1[iR+rguU]==umiArray[0]) {//other UMIs were filtered out
                uint64 iread1 = rGU1[iR+rguR];
                readInfo[iread1].cb = indCB[iCB] ;
                readInfo[iread1].umi = umiArray[0];
            };
        };
    };
    */
};
