#include "Solo.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"

void Solo::outputNumUMIperGeneCB() 
{    
    //sort by nUtot per CB
    //qsort(nUperCB,nCB,2*sizeof(uint32),funCompareNumbersReverse<uint32>); //sort by gene number

    //find 99%ile
    //uint32 nUMImax = nUperCB[2*30];//robust estimate of the max UMI
    //uint32 nUMImin = nUMImax/10;
    
    //nUMImin=0; //full output

    ofstream &geneStr=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[1],ERROR_OUT, P);
    for (uint32 ii=0; ii<Trans.nGe; ii++)
        geneStr << Trans.geID.at(ii) << '\n';
    geneStr.close();

    ofstream &cbStr=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[2],ERROR_OUT, P);
    ofstream &cbGeneMatrix=ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+pSolo.outFileNames[3],ERROR_OUT, P);

    for (auto const &cb : pSolo.cbWL)
        cbStr << convertNuclInt32toString(cb,pSolo.cbL)  <<'\n';    
    cbStr.flush();

    //header
    cbGeneMatrix <<"%%\n%\n" << Trans.nGe<< ' ' << pSolo.cbWL.size() <<' '<< nCellGeneEntries << '\n';
    
    for (uint32 icb=0; icb<nCB; icb++) {
        uint32 *rCBpp=rCBp[icb];
        for (uint32 ig=0; ig<nGperCB[icb]; ig++) {
            cbGeneMatrix << rCBpp[0]+1  <<'\t'<< indCB[icb]+1 <<'\t'<< rCBpp[1];
            if (rCBpp[1]>1) {//3 counts recorded
                cbGeneMatrix <<'\t'<< rCBpp[2] <<'\t'<< rCBpp[3];
                rCBpp += 4;
            } else {//1 recorded
                cbGeneMatrix <<"\t1\t1";
                rCBpp +=2;
            };
            cbGeneMatrix << '\n';
        };
    };
    cbGeneMatrix.flush();
       
    //*soloStatsStream << setw(50) << "Maximum number of UMIs per CB" << setw(15) << nUperCB[0] << '\n';
    //*soloStatsStream << setw(50) << "Robust maximum number of UMIs per CB" << setw(15) << nUMImax << '\n';
    //*soloStatsStream << setw(50) << "Number of CBs that passed min UMI threshold " << setw(15) << nCBout << '\n';
    //*soloStatsStream << setw(50) << "UMIs IN CELL BARCODES:\n"
    //*soloStatsStream << setw(50) << "nDetectedCellBarcodes" << setw(15) << nCB << '\n';
    
    soloStatsStream->flush();
};
