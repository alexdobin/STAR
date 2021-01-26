#include "readBarcodeLoad.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"

void loadBarcodeRead(Parameters &P, istream **readInStream, string &seq1, string &qual1)
{
    if (!P.pSolo.barcodeReadYes)
        return; //no barcode
        
    string readID;        
    getline(*readInStream[P.pSolo.barcodeRead],readID);
    getline(*readInStream[P.pSolo.barcodeRead], seq1);
    readInStream[P.pSolo.barcodeRead]->ignore(DEF_readNameSeqLengthMax,'\n');//skip to the end of 3rd ("+") line
    getline(*readInStream[P.pSolo.barcodeRead], qual1);
    
    if (seq1.size() != qual1.size()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length\n";
        errOut << readID <<'\n'<< seq1 <<'\n'<< qual1 <<'\n'<< "SOLUTION: fix your fastq file\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    
    if (seq1.size() != P.pSolo.bL) {
        if (P.pSolo.bL > 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in input read file: the total length of barcode sequence is "  << seq1.size() << " not equal to expected " <<P.pSolo.bL <<"\n"  ;
            errOut << "Read ID="<<readID<< "   Sequence="<<seq1<<"\n";
            errOut << "SOLUTION: make sure that the barcode read is the last file in --readFilesIn , and check that it has the correct formatting\n";
            errOut << "          If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if (seq1.size()<P.pSolo.cbumiL) {//barcode sequence too short - append Ns
            seq1.append(P.pSolo.cbumiL-seq1.size(), 'N');
            qual1.append(P.pSolo.cbumiL-qual1.size(), 'F');
        };
    };
    //histogram for qualities of barcode portion only
    if (P.outFilterBySJoutStage != 2) {
        g_statsAll.qualHistCalc(P.pSolo.barcodeRead, qual1.c_str()+P.pSolo.barcodeStart, P.pSolo.barcodeEnd==0 ? qual1.size() : P.pSolo.barcodeEnd-P.pSolo.barcodeStart+1);
    };
};