#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloFeature::outputResults() 
{    
    if (featureType==0) {//this only need to be done once
        //output genes
        ofstream &geneStr=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[1],ERROR_OUT, P);
        for (uint32 ii=0; ii<Trans.nGe; ii++)
            geneStr << Trans.geID[ii] <<"\t"<< (Trans.geName[ii].empty() ? Trans.geID[ii] : Trans.geName[ii]) << '\n';
        geneStr.close();
    };
    if (featureType==0 || !pSolo.featureYes[0]) {
        //output CBs
        ofstream &cbStr=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[2],ERROR_OUT, P);
        for (auto const &cb : pSolo.cbWL)
            cbStr << convertNuclInt32toString(cb,pSolo.cbL)  <<'\n';    
        cbStr.flush();
    };

    //output counting matrix
    string matrixFileName=P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[3+featureType];
    ofstream &countMatrixStream=ofstrOpen(matrixFileName,ERROR_OUT, P);

    //header
    countMatrixStream <<"%%MatrixMarket matrix coordinate integer general\n%\n";
    countMatrixStream << Trans.nGe <<' '<< pSolo.cbWL.size() <<' '<< nCellGeneEntries << '\n';
    
    for (uint32 icb=0; icb<nCB; icb++) {
        uint32 *rCBpp=rCBp[icb];
        for (uint32 ig=0; ig<nGperCB[icb]; ig++) {
            countMatrixStream << rCBpp[0]+1  <<' '<< indCB[icb]+1;
            uint32 count1[3] = {1,1,1};
            if (rCBpp[1]>1) {//3 counts recorded
                count1[0] = rCBpp[1];                
                count1[1] = rCBpp[2];
                count1[2] = rCBpp[3];                
                rCBpp += 4;
            } else {//1 recorded
                rCBpp +=2;
            };
            for (uint32 ii=0; ii<pSolo.umiDedupColumns.size(); ii++) 
                countMatrixStream <<' '<< count1[pSolo.umiDedupColumns[ii]];            
            countMatrixStream << '\n';
        };
    };
    countMatrixStream.flush();
       
    statsStream->flush();
};
