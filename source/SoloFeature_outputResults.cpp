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
            geneStr << Trans.geID[ii] <<"\t"<< Trans.geName[ii] << '\n';
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
    countMatrixStream <<"%%\n%\n" << Trans.nGe<< ' ' << pSolo.cbWL.size() <<' '<< nCellGeneEntries << '\n';
    
    for (uint32 icb=0; icb<nCB; icb++) {
        uint32 *rCBpp=rCBp[icb];
        for (uint32 ig=0; ig<nGperCB[icb]; ig++) {
            uint32 count1[3] = {1,1,1};
            count1[0] = rCBpp[1];
            if (rCBpp[1]>1) {//3 counts recorded
                count1[1]=rCBpp[2];
                count1[2]=rCBpp[3];                
                rCBpp += 4;
            } else {//1 recorded
                rCBpp +=2;
            };
            countMatrixStream << rCBpp[0]+1  <<'\t'<< indCB[icb]+1;
            for (uint32 ii=0; ii<pSolo.umiDedupColumns.size(); ii++) 
                countMatrixStream <<'\t'<< count1[pSolo.umiDedupColumns[ii]];            
            countMatrixStream << '\n';
        };
    };
    countMatrixStream.flush();
       
    //*statsStream << setw(50) << "Maximum number of UMIs per CB" << setw(15) << nUperCB[0] << '\n';
    //*statsStream << setw(50) << "Robust maximum number of UMIs per CB" << setw(15) << nUMImax << '\n';
    //*statsStream << setw(50) << "Number of CBs that passed min UMI threshold " << setw(15) << nCBout << '\n';
    //*statsStream << setw(50) << "UMIs IN CELL BARCODES:\n"
    //*statsStream << setw(50) << "nDetectedCellBarcodes" << setw(15) << nCB << '\n';
    
    statsStream->flush();
};
