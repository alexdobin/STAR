#include "SoloFeature.h"
#include "streamFuns.h"
//#include "TimeFunctions.h"
//#include "SequenceFuns.h"
//#include "Stats.h"
//#include "GlobalVariables.h"

void SoloFeature::redistributeReadsByCB()
{//redistribute reads in files by CB - each file with the approximately the same number of reads, each CB is on one file only
    
    /* SoloFeature vars that have to be setup:
     * nCB
     * readFeatSum->cbReadCount[]
    */
    
    //find boundaries for cells
    uint64 nReadRec=0;
    for ( uint32 ii=0; ii<nCB; ii++ )
        nReadRec += readFeatSum->cbReadCount[ii];
    
    uint64 nReadRecBin=nReadRec/pSolo.redistrReadsNfiles;
      
    redistrFilesCBfirst.push_back(0);
    redistrFilesCBindex.resize(nCB);
    uint64 nreads=0;
    uint32 ind=0;
    for (uint32 icb=0; icb<nCB; icb++){
        redistrFilesCBindex[icb]=ind;
        nreads += readFeatSum->cbReadCount[indCB[icb]];
        if (nreads>=nReadRecBin) {
            ind++;
            redistrFilesCBfirst.push_back(icb+1);
            redistrFilesNreads.push_back(nreads);
            nreads=0;            
        };
    };
    if (nreads>0) {
        redistrFilesCBfirst.push_back(nCB);
        redistrFilesNreads.push_back(nreads);
    };
    
    //open output files
    redistrFilesStreams.resize(redistrFilesNreads.size());
    for (uint32 ii=0; ii<redistrFilesNreads.size(); ii++)
        redistrFilesStreams[ii] = &fstrOpen(P.outFileTmp + "solo"+SoloFeatureTypes::Names[featureType]+"_redistr_"+std::to_string(ii), ERROR_OUT, P);

    //main cycle
    for (int ii=0; ii<P.runThreadN; ii++) {
        readFeatAll[ii]->streamReads->flush();
        readFeatAll[ii]->streamReads->seekg(0,ios::beg);
        
        while ( true ) {
            string line1;
            getline(*readFeatAll[ii]->streamReads,line1);
            if (line1.empty()) {
                break;
            };
            
            istringstream line1stream(line1);
            int64 cb1;            
            line1stream >> cb1 >> cb1 >> cb1;
            if (featureType==SoloFeatureTypes::SJ)
                line1stream >> cb1;
            line1stream >> cb1;
            
            *redistrFilesStreams[redistrFilesCBindex[indCBwl[cb1]]] << line1 <<'\n';
            
        };
        //TODO: delete streamReads files one by one to save disk space
    };
    
    //close files
    //for (uint32 ii=0; ii<pSolo.redistrReadsNfiles; ii++)
    //    redistrFilesStreams[ii]->flush();
};
    