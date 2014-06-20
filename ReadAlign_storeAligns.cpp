/**  ReadAlign - one read, all alignments
 *
 *
 * A longer description.
 *  
 * @see something
 */

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::storeAligns (uint iDir, uint Shift, uint Nrep, uint L, uint indStartEnd[2], uint iFrag) {//fill in alignment data

    if ( Nrep > P->seedMultimapNmax ) {// if a piece maps too many times, do not store it
        if ( Nrep < multNmin || multNmin==0 ) {multNmin=Nrep; multNminL=L;};
        return;     
    };

    nUM[ Nrep==1 ? 0:1] += Nrep;        //add numbers of U/M aligns
    nA += Nrep;
    
    uint rStart=iDir==0 ? Shift : Shift+1-L;//alignment read-start
        
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
    if (nP > P->seedPerReadNmax) {//
        ostringstream errOut;
        errOut <<"EXITING because of FATAL error: too many pieces pere read\n" ;
        errOut <<"SOLUTION: increase input parameter --seedPerReadNmax";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_RUNTIME, *P);
    };

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
