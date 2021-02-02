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
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);    
      
    typedef struct {
        uint32 feature;
        uint64 umi; //umi here are read start/end => 64bits
    } typeFeatureUMI;    
    vector<vector<typeFeatureUMI>::iterator> cbFeatUMI (nCB + 1);    
    
    #pragma omp parallel for num_threads(P.runThreadN)
    for (uint32 ired=0; ired<redistrFilesCBfirst.size()-1; ired++) {//does parallelization speed up things?
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
            
            {//collapse
                uint32 prevF = (uint32)-1;
                uint64 prevU = (uint32)-1;
                for (auto fu=cbFeatUMI[icb]; fu!=cbFeatUMI[icb]+nReadPerCB[icb]; fu++) {//cycle over all reads for this icb
                    if ( fu->feature != prevF ) {//compare to previous feature
                        prevF = fu->feature;
                        countCellGeneUMI[countCellGeneUMIindex[icb+1] + 0] = prevF;//create next feature entry
                        nGenePerCB[icb]++;
                        nUMIperCB[icb] += countCellGeneUMI[countCellGeneUMIindex[icb+1] + pSolo.umiDedup.countInd.main];
                        countCellGeneUMIindex[icb+1] = countCellGeneUMIindex[icb+1] + countMatStride;//iCB+1 accumulates the index
                    };
                    
                    if (pSolo.umiDedup.yes.NoDedup)
                        ++countCellGeneUMI[countCellGeneUMIindex[icb+1] + pSolo.umiDedup.countInd.NoDedup];
                    
                    if ( pSolo.umiDedup.yes.Exact && fu->umi != prevU ) {//same feature, new umi
                        ++countCellGeneUMI[countCellGeneUMIindex[icb+1] + pSolo.umiDedup.countInd.Exact];//collapsed UMI count
                    };
                };
            };
        };
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished reading / collapsing" <<endl;    

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
