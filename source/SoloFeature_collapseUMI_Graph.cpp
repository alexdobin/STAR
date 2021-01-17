#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include <unordered_map>
#include "SoloCommon.h"
#include <bitset>

#define def_MarkNoColor  (uint32) -1

void collapseUMIwith1MMlowHalf(uint32 *umiArr, uint32 umiArrayStride, uint32 umiMaskLow, uint32 nU0, uint32 &nU1, uint32 &nU2, uint32 &nC, vector<array<uint32,2>> &vC);
void graphDepthFirstSearch(uint32 n, vector<vector<uint32>> &nodeEdges, vector <uint32> &nodeColor);
uint32 graphNumberOfConnectedComponents(uint32 N, vector<array<uint32,2>> V, vector<uint32> &nodeColor);


uint32 SoloFeature::umiArrayCorrect_Graph(const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr)
{
    uint32 nU1 = nU0;
    uint32 nU2 = nU0;
    uint32 graphN = 0; //number of nodes
    vector<array<uint32,2>> graphConn;//node connections
    vector<uint32> graphComp;//for each node (color) - connected component number
    
    for (uint64 iu=0; iu<nU0*umiArrayStride; iu+=umiArrayStride)
        umiArr[iu+2]=def_MarkNoColor; //marks no color for graph

    qsort(umiArr, nU0, umiArrayStride*sizeof(uint32), funCompareNumbers<uint32>);
    collapseUMIwith1MMlowHalf(umiArr, umiArrayStride, pSolo.umiMaskLow, nU0, nU1, nU2, graphN, graphConn);

    //exchange low and high half of UMIs, re-sort, and look for 1MM again
    for (uint32 iu=0; iu<umiArrayStride*nU0; iu+=umiArrayStride) {
        pSolo.umiSwapHalves(umiArr[iu]);
    };
    qsort(umiArr, nU0, umiArrayStride*sizeof(uint32), funCompareNumbers<uint32>);
    collapseUMIwith1MMlowHalf(umiArr, umiArrayStride, pSolo.umiMaskLow, nU0, nU1, nU2, graphN, graphConn);

    uint32 nConnComp=graphNumberOfConnectedComponents(graphN, graphConn, graphComp);
    nU1 += nConnComp;    
    
    if (readInfoRec) {
        for (uint32 ii=0; ii<graphComp.size(); ii++) {//for non-conflicting colors, need to fill the colors correctly
            if (graphComp[ii]==(uint32)-1)
                graphComp[ii]=ii;
        };

        const uint32 bitTopMask=~(1<<31);
        vector<array<uint32,2>> umiBest(graphN,{0,0});
        unordered_map<uintUMI,uint32> umiCorrColor;
        for (uint32 iu=0; iu<umiArrayStride*nU0; iu+=umiArrayStride) {
            //switch low/high to recover original UMIs
            pSolo.umiSwapHalves(umiArr[iu]);//halves were swapped, need to return back to UMIs
            //find best UMI (highest count) for each connected component
            if (umiArr[iu+2]==def_MarkNoColor)
                continue; //UMI is not corrected
            uint32 color1=graphComp[umiArr[iu+2]];
            uint32 count1=umiArr[iu+1] & bitTopMask;
            if (umiBest[color1][0] < count1) {
                umiBest[color1][0] = count1;
                umiBest[color1][1] = umiArr[iu];
            };              
            umiCorrColor[umiArr[iu]] = color1;
        };

        for (uint32 iu=0; iu<umiArrayStride*nU0; iu+=umiArrayStride) {
            auto umi = umiArr[iu+0];
            if ( umiCorrColor.count(umi)>0 )
                umiCorr[umi]=umiBest[umiCorrColor[umi]][1];
        };
    };
    
    if (nUMIyes) {//this is not needed
        return nU1;
    } else {
        return 0;
    };
    
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
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

            //"directional" collapse from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).
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


