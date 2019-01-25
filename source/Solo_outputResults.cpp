#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void SoloFeature::outputNumUMIperGeneCB() 
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
    countMatrixStream <<"%%MatrixMarket matrix coordinate integer general\n%\n" << Trans.nGe<< ' ' << pSolo.cbWL.size() <<' '<< nCellGeneEntries << '\n';
    
    for (uint32 icb=0; icb<nCB; icb++) {
        uint32 *rCBpp=rCBp[icb];
        for (uint32 ig=0; ig<nGperCB[icb]; ig++) {
            countMatrixStream << rCBpp[0]+1  <<'\t'<< indCB[icb]+1 <<'\t'<< rCBpp[1];
            if (rCBpp[1]>1) {//3 counts recorded
                countMatrixStream <<'\t'<< rCBpp[2] <<'\t'<< rCBpp[3];
                rCBpp += 4;
            } else {//1 recorded
                countMatrixStream <<"\t1\t1";
                rCBpp +=2;
            };
            countMatrixStream << '\n';
        };
    };
    countMatrixStream.flush();
       
    //*soloStatsStream << setw(50) << "Maximum number of UMIs per CB" << setw(15) << nUperCB[0] << '\n';
    //*soloStatsStream << setw(50) << "Robust maximum number of UMIs per CB" << setw(15) << nUMImax << '\n';
    //*soloStatsStream << setw(50) << "Number of CBs that passed min UMI threshold " << setw(15) << nCBout << '\n';
    //*soloStatsStream << setw(50) << "UMIs IN CELL BARCODES:\n"
    //*soloStatsStream << setw(50) << "nDetectedCellBarcodes" << setw(15) << nCB << '\n';
    
    soloStatsStream->flush();
};
