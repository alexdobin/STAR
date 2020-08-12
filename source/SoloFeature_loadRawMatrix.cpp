#include "SoloFeature.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include "SoloFeatureTypes.h"
#include "serviceFuns.cpp"

void SoloFeature::loadRawMatrix()
{    
    //make directories
    string outputPrefix1;
    outputPrefix1 = outputPrefix + "/raw/";

    /////////////////////////////////////////////////////////////
    //output counting matrix
    string matrixFileName=outputPrefix1+pSolo.outFileNames[3];
    ifstream &matStream=ifstrOpen(matrixFileName, ERROR_OUT, "SOLUTION: check path and permission for the matrix file" + matrixFileName, P);

    //header
    while (matStream.peek() == '%') {
            string dummy;
            matStream.ignore(numeric_limits<streamsize>::max(), '\n');
    };
    
    uint32 nCB1; //number of features read from file
    uint64 nTot; //total number of entries
    matStream >> featuresNumber >> nCB1 >> nTot;
    
    /* do not need it - read the information from features file
    if (nfeat1 != featuresNumber) {
        ostringstream errOut;
        errOut <<"Exiting because of fatal INPUT FILE error: number of features in "<< matrixFileName <<": "<< nfeat1;
        errOut <<"\n   is not equal to the number of features from annotations from genome index: " << featuresNumber << "\n";
        errOut <<"SOLUTION: use the same genome index that was use in the STARsolo run that generated the matrix\n"; 
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    */
    
    //struct {uint32 gene; uint32 cb; uint32 count} matEntry;
    /*vector<array<uint32,3>> matSparse(nTot);
    for (uint64 ii=0; ii<nTot; ii++) {
        for (uint32 jj=0; jj<3; jj++)
            matStream >> matSparse[ii][jj];
    };

    //sort by cell index
    std::sort(matSparse.begin(), matSparse.end(), [](const array<uint32,3> &m1, const array<uint32,3> &m2) {
                                                        return m1[1] < m2[1] ;
                                                    });
    */
    
    countMatStride=3;
    countCellGeneUMI.resize(nTot*countMatStride);
    for (uint64 ii=0; ii<nTot; ii++) {
        matStream >> countCellGeneUMI[ii*countMatStride+0];//gene
        matStream >> countCellGeneUMI[ii*countMatStride+1];//cell
        matStream >> countCellGeneUMI[ii*countMatStride+2];//count: for now, only allow one value per cell/gene
        
        --countCellGeneUMI[ii*countMatStride+0];
        --countCellGeneUMI[ii*countMatStride+1];
    };
    
    qsort((void*) countCellGeneUMI.data(), nTot, countMatStride*sizeof(countCellGeneUMI[0]), funCompareTypeSecondFirst<uint32>);

    //count number of detected cell
    nCB=0;
    uint32 ciprev=(uint32) -1;
    for (uint32 ii=0; ii<nTot; ii++) {
        uint32 ci1=countCellGeneUMI[ii*countMatStride+1];
        if (ci1!=ciprev) {//new cell
            ciprev=ci1;
            nCB++;
        };
    };    
    
    indCB.resize(nCB);
    countCellGeneUMIindex.resize(nCB);
    nUMIperCB.resize(nCB,0);
    nGenePerCB.resize(nCB,0);
    nReadPerCB.resize(nCB,0);
    
    nCB=(uint32) -1;
    ciprev=(uint32) -1;
    for (uint32 ii=0; ii<nTot; ii++) {
        uint32 ci1 = countCellGeneUMI[ii*countMatStride+1];
        if (ci1 != ciprev) {//new cell
            ciprev = ci1;
            nCB++;
            indCB[nCB] = ci1;
            countCellGeneUMIindex[nCB] = ii*countMatStride;
        };
        nGenePerCB[nCB]++;
        nUMIperCB[nCB] += countCellGeneUMI[ii*countMatStride+2];
    };
    
    {//load barcodes
        ifstream &wlstream = ifstrOpen(outputPrefix1+pSolo.outFileNames[2], ERROR_OUT, "SOLUTION: check the path and permissions of the barcodes file", P);
        pSolo.cbWLstr.resize(nCB1);
        for (auto &cb: pSolo.cbWLstr)
            std::getline(wlstream, cb);
    };
    
    {//copy features
        
    };
    
    return;
};
