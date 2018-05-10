#include "ReadAlign.h"

void ReadAlign::waspMap() {
    if (!P.wasp.yes || trBest->varAllele.size()==0) {//no variants, vW tag will not be output
        waspType=-1;
        return;
    } else if (nTr>1) {//multimapping read
        waspType=2;
        return;
    } else if (trBest->varAllele.size()>10) {//multimapping read
        waspType=7;
        return;
    };
        
    
    waspRA->copyRead(*this);

    vector <char> vA=trBest->varAllele;

    for (const auto& a : vA) {
        if (a>3) {//read has N for the variant, drop it
            waspType=3;
            return;
        };
    };
         

    
    vector<vector<char>> vvA {{}}; //all combinations
    for (const auto& u : vA) {//cycle over vars, each time add new variant by adding 2 variants to each of the existing combinations
        (void) u; //to avoid unused warning
        vector<vector<char>> r; //temp
        for (const auto& x : vvA) {
            for (const auto& y:{1,2}) {
                r.push_back(x);
                r.back().push_back(y);
            };
        };
        vvA=move(r);
    };

    
    for (const auto& vA1 : vvA) {//cycle over all combinations

            if (vA1==vA)
                continue; //this combination was already mapped as the real read
    
            for (uint iv=0; iv<vA1.size(); ++iv) {//set all variants in this combination

                //we assume the homo-vars are already excluded
                char nt2=mapGen.Var->snp.nt[trBest->varInd.at(iv)][vA1.at(iv)]; //the other allele
                uint vr=trBest->varReadCoord.at(iv);//read coordinate

                if (trBest->Str==1) {//variant was found on the - strand alignment
                    nt2=3-nt2;
                    vr=Lread-1-vr;
                };
                waspRA->Read1[0][vr]        =nt2;
                waspRA->Read1[1][vr]        =3-nt2;
                waspRA->Read1[2][Lread-1-vr]=3-nt2;
            };

            waspRA->mapOneRead();
            waspRA->multMapSelect();
            waspRA->mappedFilter();
                    
            if (waspRA->unmapType!=-1) {
                waspType=4;
                return;
            } else if (waspRA->nTr>1) {
                waspType=5;
                return;
            } else if (waspRA->trBest->nExons!=trBest->nExons) {
                waspType=6;
                return;
            } else {
                for (uint ii=0; ii<trBest->nExons; ii++) {               
                    for (uint jj=0; jj<=2; jj++) {
                        if (trBest->exons[ii][jj]!=waspRA->trBest->exons[ii][jj]) {
                            waspType=6;
                            return;//this combination maps to a different place, return with waspType 0 (set above)
                        };
                    };
                };
            };
    }; 
    waspType=1; //all combinations resulted in the same alignment
    return;
};

void ReadAlign::copyRead(ReadAlign &r) {//copy read information only   
    Lread=r.Lread;
    readLength[0]=r.readLength[0];readLength[1]=r.readLength[1];
    readLengthOriginal[0]=r.readLengthOriginal[0];readLengthOriginal[1]=r.readLengthOriginal[1];
    readLengthPairOriginal=r.readLengthPairOriginal;
    outFilterMismatchNmaxTotal=r.outFilterMismatchNmaxTotal;
    readName=r.readName;
    
    for (uint ii=0;ii<=2;ii++)
        memcpy(Read1[ii],r.Read1[ii],Lread);//need to copy since it will be changed
    Qual1=r.Qual1;

};
