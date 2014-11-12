/**  ReadAlign - one read, all alignments
 */

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::storeAligns (uint iDir, uint Shift, uint Nrep, uint L, uint indStartEnd[2], uint iFrag) {//fill in alignment data

    #ifdef OFF_BEFORE_STORE
        #warning OFF_BEFORE_STORE
        return;
    #endif   
        
    if ( Nrep > P->seedMultimapNmax ) {// if a piece maps too many times, do not store it
        if ( Nrep < multNmin || multNmin==0 ) {multNmin=Nrep; multNminL=L;};
        return;     
    };

    nUM[ Nrep==1 ? 0:1] += Nrep;        //add numbers of U/M aligns
    nA += Nrep;
    
    uint rStart=iDir==0 ? Shift : Shift+1-L;//alignment read-start
        
  #ifdef OPTIM_STOREaligns_SIMPLE
    //find the place to insert the new entry to keep it sorted
    int iP;
    for (iP=nP-1; iP>=0; iP--) {
        if ( PC[iP][0]<=rStart ) {
            if ( (PC[iP][PC_rStart]==rStart) && PC[iP][PC_Length]<L ) continue;
            if ( (PC[iP][PC_rStart]==rStart) && PC[iP][PC_Length]==L ) return; //same alignment as before, do not store!
            break;
        };
    };

    iP=iP+1; //this is the insertion place
    for (int ii=nP-1;ii>=iP;ii--) {//move old entries to free space for the new one
        for (int jj=0;jj<PC_SIZE;jj++) PC[ii+1][jj]=PC[ii][jj];
    };

    nP++; //now nP is the new number of elements   
//     
    if (nP > P->seedPerReadNmax) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL error: too many pieces pere read\n" ;
        errOut <<"SOLUTION: increase input parameter --seedPerReadNmax";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_RUNTIME, *P);
    };    
  #else
    int iP3;
    for (iP3=nP-1; iP3>=0; iP3--) {
        if ( PC[iP3][0]<=rStart ) {
            if ( (PC[iP3][PC_rStart]==rStart) && PC[iP3][PC_Length]<L ) iP3=-1;
            if ( (PC[iP3][PC_rStart]==rStart) && PC[iP3][PC_Length]==L ) return; //same alignment as before, do not store!
            break;
        };
    };
    
    int iP2=-1;
    int nRemove=0;
    for (int iP1=0; iP1<nP; iP1++) {
        if ( (PC[iP1][PC_rStart]==rStart) && PC[iP1][PC_Length]==L ) return; //exactly the same piece
        if ( rStart >= PC[iP1][PC_rStart] ) {
            if ( rStart+L <= PC[iP1][PC_rStart]+PC[iP1][PC_Length] ) {//new piece is within the old piece
                //decide whether to keep the old one
                if ( false && (PC[iP1][PC_Nrep]>P->winAnchorMultimapNmax || Nrep<=P->winAnchorMultimapNmax)) {//old piece is not anchor, or new piece is anchor
                    PC[iP1][PC_Dir]=-1;//do not keep the old piece
                    ++nRemove;
                    if (iP2==-1) {
                        iP2=iP1; //insertion point: before iP1;
                    };                    
                } else {
                    //keep the old piece
                    if (iP2==-1) {
                        iP2=iP1+1; //insertion point: after iP1;
                    };                    
                };
            } else {//new end extends past old end
                //this is the insertion point for the new piece
                if (iP2==-1) {
                    iP2=iP1+1; //insertion point: after iP1
                };
            };
        };
        if ( rStart <= PC[iP1][PC_rStart] )  {// rStart <= PC[iP1][PC_rStart]
            if ( rStart+L >= PC[iP1][PC_rStart]+PC[iP1][PC_Length] ) {//old piece is within the new piece
                //decide whether to keep the new piece
                if ( false && (Nrep>P->winAnchorMultimapNmax || PC[iP1][PC_Nrep]<=P->winAnchorMultimapNmax) ) {//new piece is not anchor, or old piece is an anchor
                    if (nRemove>0) {//debug
                        cout << "BUG: nRemove="<<nRemove<<flush;
                        exit(-1);
                    };
                    return;//do not store the new piece
                } else {
                    //keep the new piece
                };
            } else {//old end extends past new end
                //no more overlap possible, break the cycle
                if (iP2==-1) {
                    if (iP1!=0) {//debug
                        cout << "BUG: iP1=0 iRead="<<iRead<<"   "<<readName<<flush;
                        exit(-1);
                    };  
                    iP2=iP1;
                };
                break;
            };
        };

    };

    if (nP==0) iP2=0;
    if (iP2==-1) {//debug
        cout << "BUG: iP2=-1 iRead="<<iRead<<flush;
        exit(-1);
    };   
    
    int iP=iP2;
    if (iP!=iP3+1) {
        cout << "BUG: iP!=iP3+1 iRead="<<iRead<<"   "<<readName<<flush;
        exit(-1);
    };
    
    if (nRemove==0) {//add piece
        if (nP == P->seedPerReadNmax) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL error: too many pieces pere read\n" ;
            errOut <<"SOLUTION: increase input parameter --seedPerReadNmax";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_RUNTIME, *P);
        };        
        for (int ii=nP-1;ii>=iP;ii--) {//move old entries to free space for the new one
            for (int jj=0;jj<PC_SIZE;jj++) PC[ii+1][jj]=PC[ii][jj];
        };
        ++nP;
    } else if (nRemove>1) {//replace in place
        int iP3=iP2+1;
        for (int ii=iP2+1; ii<nP; ii++) {//move old entries to free space for the new one
            if ( PC[ii][PC_Dir] != -1) {//move this record to the new place
                for (int jj=0;jj<PC_SIZE;jj++) PC[iP3][jj]=PC[ii][jj];//TODO use memcopy
                ++iP3;
            };
        };        
        nP-=nRemove-1;
    };
  #endif
    
    //store new piece
    PC[iP][PC_rStart]=rStart;  //alignment read-start     
    PC[iP][PC_Length]=L;       //alignment length
    PC[iP][PC_Dir]    = iDir; //direction        
    PC[iP][PC_Nrep]   = Nrep; //repeat number - for both strands        
    PC[iP][PC_SAstart]= indStartEnd[0]; //SA index 1
    PC[iP][PC_SAend]  = indStartEnd[1]; //SA index 2       
    PC[iP][PC_iFrag]  = iFrag;
   
    //choose "best" alignment

    if (L<storedLmin) L=storedLmin;

    if (Nrep==1) {
        if (L>uniqLmax) {
            uniqLmax=L;
            uniqLmaxInd=nP-1;
        };
    } else {
        if ( Nrep < multNmin || multNmin==0 ) {multNmin=Nrep; multNminL=L;};
        if ( L > multLmax ) {multLmax=L;multLmaxN=Nrep;};
        if ( Nrep > multNmax ) {multNmax=Nrep; multNmaxL=L;};
    };
};
