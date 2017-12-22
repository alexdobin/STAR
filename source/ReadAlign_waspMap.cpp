#include "ReadAlign.h"

void ReadAlign::waspMap() {
    waspType=-1;
    if (!P.wasp.yes || nTr>1 || trBest->varAllele.size()==0)
        return;
    
    waspRA->copyRead(*this);
    waspType=0;
    
    if (trBest->varAllele.size()==1) {
        
        int iv=0;

        //we assume the homo-vars are already excluded
        char nt2=mapGen.Var->snp.nt[trBest->varInd.at(iv)][3-trBest->varAllele.at(iv)]; //the other allele
        uint vr=trBest->varReadCoord.at(iv);//read coordinate
        
        if (trBest->Str==1) {//variant was found on the - strand alignment
            nt2=3-nt2;
            vr=Lread-1-vr;
        };
        waspRA->Read1[0][vr]        =nt2;
        waspRA->Read1[1][vr]        =3-nt2;
        waspRA->Read1[2][Lread-1-vr]=3-nt2;
                
        waspRA->mapOneRead();
        
        if (waspRA->unmapType==-1 && waspRA->nTr==1 && waspRA->trBest->nExons==trBest->nExons) {
            for (uint ii=0; ii<trBest->nExons; ii++) {
//                 uint r=
                
                for (uint jj=0; jj<=2; jj++) {
                    if (trBest->exons[ii][jj]!=waspRA->trBest->exons[ii][jj])
                        return;
                };
            };
            waspType=1;
        };
    } else {
        waspType=2;
    };
    
    return;
};

void ReadAlign::copyRead(ReadAlign &r) {//copy read information only   
    Lread=r.Lread;
    readLength[0]=r.readLength[0];readLength[1]=r.readLength[1];
    readLengthOriginal[0]=r.readLengthOriginal[0];readLengthOriginal[1]=r.readLengthOriginal[1];
    readLengthPairOriginal=r.readLengthPairOriginal;
    outFilterMismatchNmaxTotal=r.outFilterMismatchNmaxTotal;

    for (uint ii=0;ii<=2;ii++)
        memcpy(Read1[ii],r.Read1[ii],Lread);//need to copy since it will be changed
    Qual1=r.Qual1;

};