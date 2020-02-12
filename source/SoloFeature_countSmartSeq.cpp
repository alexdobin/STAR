#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include "soloInputFeatureUMI.h"


void SoloFeature::countSmartSeq()
{    
    time_t rawTime;
    
    //nCB=pSolo.cbWLsize; //all cells are recorded
    //nCB, indCB are recorded in sumThreads

    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatSum->addStats(*readFeatAll[ii]);
    };
    
    redistributeReadsByCB();
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ... Finished redistribution of reads from Solo read files"<<endl;   
    
    //////////////////////////////////////////////////////
    nReadPerCB.resize(nCB);
    cbFeatureUMImap.resize(nCB);
    
    //rguStride=3;//uint32 for gene, uint64 for UMI
    
    typedef struct {
        uint32 feature;
        uint32 count[2];
    } typeFeatureCount;   
    
    vector<vector<typeFeatureCount>> vCellFeatureCount(nCB);
    
    typedef struct {
        uint32 feature;
        uint64 umi;
    } typeFeatureUMI;    
    vector<vector<typeFeatureUMI>::iterator> cbFeatUMI (nCB + 1);    
    
    #pragma omp parallel for num_threads(P.runThreadN)
    for (uint32 ired=0; ired<redistrFilesCBfirst.size()-1; ired++) {//TODO: parallelize, each ired is independent here!
        vector<typeFeatureUMI> vFeatureUMI (redistrFilesNreads[ired]);
        
        uint32 iCB1=redistrFilesCBfirst[ired];
        uint32 iCB2=redistrFilesCBfirst[ired+1];
        
        //allocate arrays    
        cbFeatUMI[iCB1]=vFeatureUMI.begin();
        for (uint32 icb=iCB1+1; icb<iCB2; icb++) {
            cbFeatUMI[icb] = cbFeatUMI[icb-1] + readFeatSum->cbReadCount[indCB[icb-1]];
        };        
        
        //input records
        redistrFilesStreams[ired]->flush();
        redistrFilesStreams[ired]->seekg(0,ios::beg);
        
        uint32 feature;
        uint64 umi, iread;
        int32 cbmatch;
        int64 cb;
        vector<uint32> trIdDist; //not used
        bool readInfoYes=false;
        while (soloInputFeatureUMI(redistrFilesStreams[ired], featureType, readInfoYes, P.sjAll, iread, cbmatch, feature, umi, trIdDist)) {//cycle over file records
            
            *redistrFilesStreams[ired] >> cb;
            
            if (feature+1 == 0) 
                continue;

            uint32 icb=indCBwl[cb];
            *( cbFeatUMI[icb] + nReadPerCB[icb] )={feature,umi};
            nReadPerCB[icb]++;
        };
        
        //collapse UMI, simple
        for (uint32 icb=iCB1; icb<iCB2; icb++) {
            if (nReadPerCB[icb]==0) {
                continue;
            };
            
            sort(cbFeatUMI[icb], cbFeatUMI[icb]+nReadPerCB[icb], 
                    [](const typeFeatureUMI &a, const typeFeatureUMI &b) 
                      {return (a.feature < b.feature) || ( a.feature == b.feature && a.umi < b.umi); });
                  
            vCellFeatureCount[icb].reserve(8192);
            vCellFeatureCount[icb].push_back({cbFeatUMI[icb]->feature, {1,1}});
            for (auto fu=cbFeatUMI[icb]+1; fu!=cbFeatUMI[icb]+nReadPerCB[icb]; fu++) {//cycle over all reads for this icb
                if ( fu->feature != (fu-1)->feature ) {//compare to previous feature
                    vCellFeatureCount[icb].push_back({fu->feature, {1,0}});//create next feature entry
                } else if ( fu->umi != (fu-1)->umi ) {//same feature, new umi
                    vCellFeatureCount[icb].back().count[0]++;//add 1
                };
                vCellFeatureCount[icb].back().count[1]++;       
            };
        };
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading / collapsing" <<endl;    
    
    //convert to countCellGeneUMI. TODO - this is not necessary, recode output to use vCellFeatureCount
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
      
    countMatStride=3;
    
    uint64 ccgN=0;//total number of entries in the countCellGeneUMI
    for (auto &cbf :  vCellFeatureCount) {
        ccgN+=cbf.size();
    };
    countCellGeneUMI.resize(ccgN*countMatStride);
    countCellGeneUMIindex.resize(nCB+1);
    countCellGeneUMIindex[0]=0;
    
    for (uint32 icb=0; icb<nCB; icb++) {//copy vCellFeatureCount
                
        //this assumes vCellFeatureCount and countCellGeneUMI have the same stride
        memcpy( (void*) (countCellGeneUMI.data() + countCellGeneUMIindex[icb]), (void*) vCellFeatureCount[icb].data(),  vCellFeatureCount[icb].size()*countMatStride*sizeof(countCellGeneUMI[0]));
        countCellGeneUMIindex[icb+1] = countCellGeneUMIindex[icb] + vCellFeatureCount[icb].size()*countMatStride;
        
        nGenePerCB[icb]=vCellFeatureCount[icb].size();
        for (auto &cbf :  vCellFeatureCount[icb]) {
            nUMIperCB[icb]+=cbf.count[pSolo.umiDedupColumns[0]];
        };
    };

    // sum stats
    uint32 nReadPerCBmax=0;
    for (uint32 icb=0; icb<nCB; icb++) {
        
        nReadPerCBmax=max(nReadPerCBmax,nReadPerCB[icb]);
        readFeatSum->stats.V[readFeatSum->stats.nMatch] += nReadPerCB[icb];
                
        readFeatSum->stats.V[readFeatSum->stats.nUMIs] += nUMIperCB[icb];
        if (nGenePerCB[icb]>0)
            ++readFeatSum->stats.V[readFeatSum->stats.nCellBarcodes];        
    };
    
    readFeatSum->stats.V[readFeatSum->stats.nExactMatch]=readFeatSum->stats.V[readFeatSum->stats.nMatch];            
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished SmartSeq counting" <<endl;
};
