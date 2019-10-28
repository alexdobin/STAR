#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "SoloCommon.h"
#include <unordered_map>
#include <bitset>
#include "serviceFuns.cpp"

void SoloFeature::quantTranscript()
{//velocyto counting gets info from Gene counting
    time_t rawTime;

    if (pSolo.clusterCBfile=="-")
        return;//no transcript quantification w/o cluster file
    
    std::unordered_map<uint32,uint32> clusterCBind; //cluster index of detected CB
    {//load cluster information
        ifstream &clusterStream = ifstrOpen(pSolo.clusterCBfile, ERROR_OUT, "SOLUTION: check the path and permissions of the cluster CB file: " + pSolo.clusterCBfile, P);
        //std::set<uint32> clusterInd; //all cluster indexes
        string seq1;
        while (clusterStream >> seq1) {
            uint64 cb1;
            if (convertNuclStrToInt64(seq1,cb1)) {//convert to 2-bit format
                uint64 cb1=binarySearchExact<uint64>(cb1,pSolo.cbWL.data(),pSolo.cbWLsize);
                if (cb1+1!=0) {
                    //uint32 icl;
                    clusterStream >> clusterCBind[cb1];
                } else {
                    P.inOut->logMain << "WARNING: cluster CB sequence not present in whitelist and is ignored: " << seq1 <<endl;
                };
            } else {
                P.inOut->logMain << "WARNING: cluster CB sequence contains non-ACGT base and is ignored: " << seq1 <<endl;
            };        
        };
    };
    
    typedef struct {
        uint32 tr;
        double d;
    } transcriptDistProbStruct;
    
    map<uint32, unordered_map<typeUMI,vector<transcriptDistProbStruct>>> mapTrDist;
    
    //////////// input records
    for (int iThread=0; iThread<P.runThreadN; iThread++) {//TODO: this can be parallelized
        fstream *streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0,ios::beg);    
        
        string line1;
        while (std::getline(*streamReads, line1)) {//until the end of file
//             stringstream lineStream(line1);
//             uint32 cb, umi, nTr;
//             lineStream >> cb >> umi >> nTr;
// 
//             if (clusterCBind.count(cb)==0) //this cb is not in the clusters
//                 continue;
//                    
//             vector<transcriptDistProbStruct> tD(nTr);
//             for (auto & tt: tD) {
//                 lineStream >> tt.tr >> tt.d;
//             }; 
//             
//             if (mapTrDist[cb].count(umi)==0) {//1st entry for this umi
//                 mapTrDist[cb][umi]=tD;
//                 continue;
//             };
// 
//             uint32 inew=0;
//             vector<trTypeStruct> tD1;
//             tD1.reserve(mapTrDist[cb][umi].size());
// 
//             for (uint32 iold=0; iold<mapTrDist[cb][umi].size(); iold++) {//intersection of old with new
//                 while (inew < tD.size() && mapTrDist[cb][umi][iold].tr>tD[inew].tr) //move through the sorted lists
//                     ++inew;
// 
//                 if (inew == tD.size() ) //end of tD reached
//                     break;
// 
//                 if (mapTrDist[cb][umi][iold].tr == tD[inew].tr) {//
//                     tD1.push_back({tD[inew].tr, (uint8)(mapTrDist[cb][umi][iold].type | tD[inew].type)});
//                 };
//             };
//             mapTrDist[cb][umi]=tD1;//replace with intersection
        };        
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Transcript3p counting: finished input" <<endl;            
 
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Transcript3p counting: finished transcript quantification" <<endl;
};
