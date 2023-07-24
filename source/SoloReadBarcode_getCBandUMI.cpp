#include "SoloReadBarcode.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#include <chrono>
#include <thread>

void SoloReadBarcode::matchCBtoWL(string &cbSeq1, string &cbQual1, vector<uint64> &cbWL, int32 &cbMatch1, vector<uint64> &cbMatchInd1, string &cbMatchString1)
{
    cbMatch1=-1;
    cbMatchString1="";
    cbMatchInd1.clear();
    //convert CB and check for Ns
    uint64 cbB1;
    int64 posN=convertNuclStrToInt64(cbSeq1,cbB1);

    if (!pSolo.cbWLyes) {//no whitelist - no search
        if (posN!=-1) {//Ns are present, discard this read
            cbMatch1=-2;
            //stats.V[stats.nNinBarcode]++;
        } else {//no Ns
            //cbI=(int64) cbB1;
            cbMatchInd1.push_back(cbB1);//all possible barcodes are accepted. This will overflow if CB is longer than 31b
            cbMatchString1 = to_string(cbB1);
            cbMatch1=0;
        };
        return;
    };

    if (posN==-2) {//>2 Ns, might already be filtered by Illumina
        cbMatch1=-2;
        //stats.V[stats.nNinBarcode]++;
        return;
    } else if (posN==-1) {//no Ns, count only for featureType==gene
        int64 cbI=binarySearchExact<uint64>(cbB1,cbWL.data(),cbWL.size());
        if (cbI>=0) {//exact match
            cbMatchInd1.push_back((uint64) cbI);
            cbMatchString1 = to_string(cbMatchInd1[0]);
            cbMatch1=0;
            return;
        };
    };
    
    if (!pSolo.CBmatchWL.mm1) //only exact matches allowed
        return;

    cbMatch1=0;
    if (posN>=0) {//one N
        uint32 posNshift=2*(cbSeq1.size()-1-posN);//shift bits for posN
        bool matched = false;
        for (uint32 jj=0; jj<4; jj++) {
            uint64 cbB11=cbB1^(jj<<posNshift);
            int64 cbI1=binarySearchExact<uint64>(cbB11,cbWL.data(),cbWL.size());
            if (cbI1>=0) {//found match
                if (!pSolo.CBmatchWL.mm1_multi_Nbase && matched) {
                    cbMatchInd1.clear();
                    cbMatch1=-3;
                    break; //this is 2nd match, not allowed for N-bases
                };
                matched = true;
                //output all
                cbMatchInd1.push_back(cbI1);
                ++cbMatch1;
                cbMatchString1 += ' ' +to_string(cbI1) + ' ' + cbQual1[posN];
            };
        };
    } else {//look for 1MM; posN==-1, no Ns
        for (uint32 ii=0; ii<cbSeq1.size(); ii++) {
            for (uint32 jj=1; jj<4; jj++) {
                int64 cbI1=binarySearchExact<uint64>(cbB1^(jj<<(ii*2)),cbWL.data(),cbWL.size());
                if (cbI1>=0) {//found match
                    //output all
                    cbMatchInd1.push_back(cbI1);
                    ++cbMatch1;
                    cbMatchString1 += ' ' +to_string(cbI1) + ' ' + cbQual1.at(cbSeq1.size()-1-ii);
                };
            };
        };
    };
    if (cbMatch1==0) {//no matches
        //stats.V[stats.noNoWLmatch]++;
        cbMatch1=-1;
    } else if (cbMatch1==1) {//1 match, no need to record the quality
        cbMatchString1 = to_string(cbMatchInd1[0]);
    } else if (!pSolo.CBmatchWL.mm1_multi) {//>1 matches, but this is not allowed
        cbMatch1=-3;
        cbMatchInd1.clear();
        cbMatchString1="";
    };// else cbMatch contains number of matches, and cbMatchString has CBs and qualities
};

void SoloReadBarcode::addStats(const int32 cbMatch1)
{
    if (!pSolo.cbWLyes) //no stats if no WL
        return;
       
    switch (cbMatch1) {
        case 0://exact matches
            cbReadCountExact[cbMatchInd[0]]++;//note that this simply counts reads per exact CB, no checks of genes or UMIs
            stats.V[stats.yesWLmatchExact]++;
            break;
        case 1: //one WL match counted here, but they may still get rejected in SoloReadFeature_inputRecords.cpp
            stats.V[stats.yesOneWLmatchWithMM]++;
            break;
        default: //multiple WL matches are counted here, but they may still get rejected in SoloReadFeature_inputRecords.cpp
            stats.V[stats.yesMultWLmatchWithMM]++;
            break;            
        case -1 :
            stats.V[stats.noNoWLmatch]++;
            break;
        case -2 :
            stats.V[stats.noNinCB]++;
            break;
        case -3 :
            stats.V[stats.noTooManyWLmatches]++;
            break;
        case -11 :            
            stats.V[stats.noNoCB]++;//CB sequence cannot be extracted
            break;
        case -12 :
            stats.V[stats.noTooManyMM]++;            
            break;
        case -23:
            stats.V[stats.noNinUMI]++;//UMIs are not allowed to have Ns
            break;
        case -24:
            stats.V[stats.noUMIhomopolymer]++;            
    };
};

///////////////////////////////////////////////////////////////////////////////////////
bool SoloReadBarcode::convertCheckUMI()
{//check UMIs, return if bad UMIs
    if (convertNuclStrToInt64(umiSeq,umiB)!=-1) {//convert and check for Ns
        umiCheck=-23;
        return false;
    };
    if (umiB==homoPolymer[0] || umiB==homoPolymer[1] || umiB==homoPolymer[2] || umiB==homoPolymer[3]) {
        umiCheck=-24;        
        return false;
    };
    return true;
};

//////////////////////////////////////////////////////////////////////////////////////
void SoloReadBarcode::getCBandUMI(char **readSeq, char **readQual, uint64 *readLen, const string &readNameExtraIn, const uint32 &readFilesIndex, const char *readName)
{
    if (pSolo.type==0)
        return;

///////////////////////////SmartSeq  
    if (pSolo.type==pSolo.SoloTypes::SmartSeq) {        
        cbSeq=cbQual=cbSeqCorrected=""; //TODO make cbSeq=file label
        cbMatch=0;
        cbMatchInd={readFilesIndex};
        cbMatchString=to_string(cbMatchInd[0]);
        addStats(cbMatch);
        return;
    };    
    
    cbMatch=-1;
    cbMatchString="";
    cbMatchInd.clear();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////    
    ////////// bSeq and bQual
    if (P.readFilesTypeN != 10) {//not SAM: barcode seq/qual are at the beginning of readNameExtra
        bSeq=std::string(readSeq[pSolo.barcodeRead], readLen[pSolo.barcodeRead]);
        bQual=std::string(readQual[pSolo.barcodeRead], readLen[pSolo.barcodeRead]);
        
    } else {//SAM: barcode seq/qual is collected from SAM tags
        uint32 pos1 = readNameExtraIn.find_first_not_of(" \t"); //to remove whitespace
        string readNameExtraT;
        if (pos1==0) {
            readNameExtraT = '\t' + readNameExtraIn;//\t needs to be in front of readNameExtra for efficient search
        } else {
            readNameExtraT = readNameExtraIn;
            readNameExtraT[pos1] = '\t';
        };
        
        bSeq = {};
        bStrings.clear();
        for (auto &tag: pSolo.samAtrrBarcodeSeq) {
            size_t pos1 = readNameExtraT.find(tag); //find tag
            if ( pos1 == std::string::npos ) {
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR in input read file: could not find barcode sequence SAM attribute "  << tag << " in read " <<readName <<"\n" ;
                errOut << "with SAM attributes: "<< readNameExtraT <<"\n";
                errOut << "SOLUTION: make sure that all reads in the input SAM/BAM have all attributes from --soloInputSAMattrBarcodeSeq\n";
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);                
            };
            pos1 +=6; //skip 6 chars, e.g. \tCB:Z:
            size_t pos2 = readNameExtraT.find('\t', pos1); //find next \t
            bStrings.push_back(readNameExtraT.substr(pos1,pos2-pos1));
            bSeq += bStrings.back();
            
        };

        bQual = {};        
        if (pSolo.samAtrrBarcodeQual.size()==0) {//if quality tags are not supplied
            bQual.resize(bSeq.size(), 'H');
        } else {
            for (auto &tag: pSolo.samAtrrBarcodeQual) {
                size_t pos1 = readNameExtraT.find(tag); //find tag, and skip 6 chars, e.g. \tCB:Z:
                if ( pos1 == std::string::npos ) {
                    ostringstream errOut;
                    errOut << "EXITING because of FATAL ERROR in input read file: could not find barcode qualities SAM attribute "  << tag << " in read " <<readName <<"\n" ;
                    errOut << "with SAM attributes: "<< readNameExtraT <<"\n";
                    errOut << "SOLUTION: make sure that all reads in the input SAM/BAM have all attributes from --soloInputSAMattrBarcodeQual\n";
                    exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);                
                };      
                pos1 += 6; //skip 6 chars, e.g. \tCB:Z:
                size_t pos2 = readNameExtraT.find('\t', pos1); //find next \t
                bQual += readNameExtraT.substr(pos1,pos2-pos1);
            };
        };
                
        if (bQual.size() != bSeq.size()) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in input read file: the total length of barcode qualities is "  << bQual.size() << " not equal to the sequence length " << bSeq.size() <<"\n"  ;
            errOut << "Read ID="<< readName <<" ;  Qualities="<< bQual <<" ;  Sequence="<< bSeq << " ;  Read SAM attributes: "<< readNameExtraT <<"\n";
            errOut << "SOLUTION: make sure correct attributes are listed in --soloInputSAMattrBarcodeQual\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };            
    };
    
    if ( bSeq.size() != P.pSolo.bL ) {//check barcodeRead length. bSeq.size == bQual.size here, this should have been checked before
        if (P.pSolo.bL > 0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in input read file: the total length of barcode sequence is "  << bSeq.size() << " not equal to expected " << P.pSolo.bL <<"\n" <<
                      "Read ID=" << readName <<" ;  Sequence=" << bSeq <<'\n' <<
                      "SOLUTION: check the formatting of input read files.\n" <<
                      "If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength\n" <<
                      "To avoid checking of barcode read length, specify --soloBarcodeReadLength 0" ;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        } else if (bSeq.size()<P.pSolo.cbumiL) {//barcode sequence too short - append Ns, will work only if P.pSolo.cbumiL>0
            bSeq.append(P.pSolo.cbumiL-bSeq.size(), 'N');
            bQual.append(P.pSolo.cbumiL-bQual.size(), 'H');
        };
    };

    if ( pSolo.type != pSolo.SoloTypes::CB_UMI_Simple ) {
        for (uint64 ix=0; ix<( P.pSolo.bL>0 ? P.pSolo.bL : bQual.size() ); ix++) {//bL==0 use the whole barcode read for quality scores histogram
            qualHist[(uint8)bQual[ix]]++;
        };
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// cbSeq and umiSeq, and match CB to WL, different CB_UMI types
    ///////////////////////////CB_UMI_Simple a.k.a Droplet
    if ( pSolo.type==pSolo.SoloTypes::CB_UMI_Simple ) {

        if (pSolo.CBtype.type==1) {//sequence cb
            cbSeq=bSeq.substr(pSolo.cbS-1,pSolo.cbL);
            umiSeq=bSeq.substr(pSolo.umiS-1,pSolo.umiL);
            cbQual=bQual.substr(pSolo.cbS-1,pSolo.cbL);
            umiQual=bQual.substr(pSolo.umiS-1,pSolo.umiL);


            for (uint64 ix=0; ix<cbQual.size(); ix++) {
                qualHist[(uint8)cbQual[ix]]++;
            };
            for (uint64 ix=0; ix<umiQual.size(); ix++) {
                qualHist[(uint8)umiQual[ix]]++;
            };               
            
            matchCBtoWL(cbSeq, cbQual, pSolo.cbWL, cbMatch, cbMatchInd, cbMatchString);
        } else if (pSolo.CBtype.type==2) {//string cb
            /* this seg-faults
            while (pSolo.CBtype.strMap.count(cbSeq)==0) {
                if (pSolo.CBtype.strMtx->try_lock()) {
                    pSolo.CBtype.strMap[cbSeq] = pSolo.CBtype.strMap.size();
                    pSolo.CBtype.strMtx->unlock();
                };
                    //std::this_thread::sleep_for(std::chrono::nanoseconds(10));
            };
            */
            cbSeq = bStrings[0];
            umiSeq = bStrings[1];
            pSolo.CBtype.strMtx->lock();
            auto cbi = pSolo.CBtype.strMap.find(cbSeq);
            uint32 cb1;
            if (cbi == pSolo.CBtype.strMap.end()) {
                cb1 = pSolo.CBtype.strMap.size();
                pSolo.CBtype.strMap[cbSeq] = cb1;
            } else {
                cb1 = cbi->second;
            };
            pSolo.CBtype.strMtx->unlock();
            
            cbMatchInd.push_back(cb1);//all possible barcodes are accepted. This will overflow if CB is longer than 31b
            cbMatchString = std::to_string(cb1);
            cbMatch=0;
        };

        if (!convertCheckUMI()) {//UMI conversion
            #ifdef MATCH_CellRanger
            // this is what CR does - will not use it to avoid differences with / 2.7.6a
            // this affects <10 reads raw, only 1 UMI count filtered
            if (cbMatch==0)
                cbReadCountExact[cbMatchInd[0]]++; //still need to count it as exact before return, even if UMI is not good
            #endif
            cbMatch=umiCheck;
            cbMatchString="";
            cbMatchInd.clear();
            addStats(cbMatch);
            return;
        };

    ///////////////////////////CB_samTagOut
    } else if ( pSolo.type==pSolo.SoloTypes::CB_samTagOut ) {//similar to CB_UMI_Simple, but no UMI, and define cbSeqCorrected
        cbSeq=bSeq.substr(pSolo.cbS-1,pSolo.cbL);
        umiSeq=bSeq.substr(pSolo.umiS-1,pSolo.umiL);
        cbQual=bQual.substr(pSolo.cbS-1,pSolo.cbL);
        umiQual=bQual.substr(pSolo.umiS-1,pSolo.umiL);

        matchCBtoWL(cbSeq, cbQual, pSolo.cbWL, cbMatch, cbMatchInd, cbMatchString);

        if ( cbMatch==0 || cbMatch==1 ) {
            if (pSolo.cbWLyes) {
        	    cbSeqCorrected = pSolo.cbWLstr[cbMatchInd[0]];
            } else {
                cbSeqCorrected = cbSeq; //no WL - no correction
            };
        } else {
        	cbSeqCorrected="-";
        };
        
    ///////////////////////////CB_UMI_Complex
    } else if ( pSolo.type==pSolo.SoloTypes::CB_UMI_Complex ) {
        
        cbSeq="";
        cbQual="";
        umiSeq="";
        umiQual="";
        
        uint32 adapterStart=0;
        if (pSolo.adapterYes) {
            if (localAlignHammingDist(bSeq, pSolo.adapterSeq, adapterStart) > pSolo.adapterMismatchesNmax) {
                stats.V[stats.noNoAdapter]++;
                cbMatch=-21;
                return; //no adapter found
            };
        };

        if (!pSolo.umiV.extractBarcode(bSeq, bQual, adapterStart, umiSeq, umiQual)) {
            stats.V[stats.noNoUMI]++;
            cbMatch=-22;
            return;
        };

        if ( pSolo.umiL == 0 )
            pSolo.umiL = umiSeq.size();

        bool cbMatchGood=true;
        if (!convertCheckUMI()) {
            cbMatchGood = false;//CB matching will not be done, just extract the sequences
            cbMatch = umiCheck;
        };

        cbMatchInd={0};
        for (auto &cb : pSolo.cbV) {//cycle over multiple barcodes
            
            string cbSeq1, cbQual1;
            if ( !cb.extractBarcode(bSeq, bQual, adapterStart, cbSeq1, cbQual1) 
                 || cbSeq1.size() < cb.minLen || cbSeq1.size() >= cb.wl.size() || cb.wl[cbSeq1.size()].size()==0 ) {
                //barcode cannot be extracted
                if (cbMatchGood) {
                    cbMatch=-11;
                    cbMatchGood=false;
                };
            };
            cbSeq  += cbSeq1 + "_";
            cbQual += cbQual1 + "_";
            
            if (!cbMatchGood)
                continue; //continue - to be able to record full cbSeq, cbQual, but no need to match to the WL
            
            uint64 cbLen1 = cbSeq1.size();
            if ( pSolo.CBmatchWL.EditDist_2 ) {
                cbMatch=0; //0: no MM, 1: >=1MM. 2: multi-match: not allowed for this option
                uintCB cbB1;
                int64 posN=convertNuclStrToInt64(cbSeq1,cbB1);
                if (posN!=-1) {//Ns in barcode: no good match
                    cbMatch = -2;
                    cbMatchGood = false;
                } else {
                    int64 cbI=binarySearchExact<uint64>(cbB1,cb.wl[cbSeq1.size()].data(),cb.wl[cbLen1].size());
                    if (cbI>=0) {//exact match
                        cbMatchInd[0] += cb.wlFactor*(cbI+cb.wlAdd[cbLen1]);
                    } else {//no exact match
                        cbI=binarySearchExact<uint64>(cbB1,cb.wlEd[cbLen1].data(),cb.wlEd[cbLen1].size());
                        if (cbI>=0) {//find match in the edited list
                            cbMatch = 1; //>=1MM
                            cbI = cb.wlEdInd[cbLen1][cbI];
                            cbMatchInd[0] += cb.wlFactor*(cbI+cb.wlAdd[cbLen1]);
                        } else {
                            cbMatch = -1;
                            cbMatchGood = false;
                        };
                    };
                };
            } else {// Exact or 1MM
                int32 cbMatch1;
                vector<uint64> cbMatchInd1;
                matchCBtoWL(cbSeq1, cbQual1, cb.wl[cbLen1], cbMatch1, cbMatchInd1, cbMatchString); //cbMatchString is not used for now, multiple matches are not allowed
                if (cbMatch1<0) {//no match
                    cbMatchGood=false;
                    cbMatch = cbMatch1;
                } else if (cbMatch1>0 && cbMatch>0) {//this barcode has >1 1MM match, or previous barcode had a mismatch
                    cbMatchGood=false;
                    cbMatch = -12; //marks mismatches in multiple barcodes
                } else {//good match
                    cbMatchInd[0] += cb.wlFactor*(cbMatchInd1[0]+cb.wlAdd[cbLen1]);
                    cbMatch=max(cbMatch,cbMatch1);//1 wins over 0
                };
            };
        };
        cbSeq.pop_back();//remove last "_" from file
        cbQual.pop_back();
        
        if (cbMatchGood) {
            cbMatchString=to_string(cbMatchInd[0]);
        };
    };
    
    addStats(cbMatch);
};
