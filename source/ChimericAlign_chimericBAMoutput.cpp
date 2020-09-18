#include "ChimericAlign.h"
#include "ReadAlign.h"
#include "BAMfunctions.h"

#include <vector>

void ChimericAlign::chimericBAMoutput(Transcript *al1, Transcript *al2, ReadAlign *RA, const uint iTr, const uint chimN, const bool isBestChimAlign, const Parameters& P)
{
    vector<Transcript*> trChim(2);
    trChim[0] = al1;
    trChim[1] = al2;

    int chimRepresent=-999, chimType=0;
    if (trChim[0]->exons[0][EX_iFrag]!=trChim[0]->exons[trChim[0]->nExons-1][EX_iFrag]) {//tr0 has both mates
        chimRepresent = 0;
        chimType = 1;
    } else if (trChim[1]->exons[0][EX_iFrag]!=trChim[1]->exons[trChim[1]->nExons-1][EX_iFrag]) {//tr1 has both mates
        chimRepresent = 1;
        chimType = 1;
    } else if (trChim[0]->exons[0][EX_iFrag]!=trChim[1]->exons[0][EX_iFrag]) {//tr0 and tr1 are single different mates
        chimRepresent = -1;
        chimType = 2;
    } else  {//two chimeric segments are on the same mate - this can only happen for single-end reads
        chimRepresent = (trChim[0]->maxScore > trChim[1]->maxScore) ? 0 : 1;
        chimType = 3;
    };

    int alignType, bamN=0, bamIsuppl=-1, bamIrepr=-1;
    uint bamBytesTotal=0;//estimate of the total size of all bam records, for output buffering
    uint mateChr,mateStart;
    uint8_t mateStrand;
    for (uint itr=0;itr<trChim.size();itr++) {//generate bam for all chimeric pieces
        trChim[itr]->primaryFlag=isBestChimAlign;
        if (chimType==2) {//PE, encompassing
            mateChr=trChim[1-itr]->Chr;
            mateStart=trChim[1-itr]->exons[0][EX_G];
            mateStrand=(uint8_t) (trChim[1-itr]->Str!=trChim[1-itr]->exons[0][EX_iFrag]);
            alignType=-10;
        } else {//spanning chimeric alignment, could be PE or SE
            mateChr=-1;mateStart=-1;mateStrand=0;//no need fot mate info unless this is the supplementary alignment
            if (chimRepresent==(int)itr) {
                alignType=-10; //this is representative part of chimeric alignment, record is as normal; if encompassing chimeric junction, both are recorded as normal
                bamIrepr=bamN;
                if (trChim[itr]->exons[0][EX_iFrag]!=trChim[1-itr]->exons[0][EX_iFrag]) {//the next mate is chimerically split
                    ++bamIrepr;
                };
//TODO check flags of SE split read
//                if (chimType==3) {
//                    bamIrepr=bamN;
//                } else if (trChim[itr]->exons[0][EX_iFrag]==trChim[1-itr]->exons[0][EX_iFrag]) {
//
//                };
//                bamIrepr=( (itr%2)==(trChim[itr]->Str) && chimType!=3) ? bamN+1 : bamN;//this is the mate that is chimerically split
            } else {//"supplementary" chimeric segment
                alignType=P.pCh.out.bamHardClip ? ( ( itr%2==trChim[itr]->Str ) ? -12 : -11) : -13 ; //right:left chimeric junction
                bamIsuppl=bamN;
                if (chimType==1) {//PE alignment, need mate info for the suppl
                    uint iex=0;
                    for (;iex<trChim[chimRepresent]->nExons-1;iex++) {
                        if (trChim[chimRepresent]->exons[iex][EX_iFrag]!=trChim[itr]->exons[0][EX_iFrag]) {
                            break;
                        };
                    };
                    mateChr=trChim[chimRepresent]->Chr;
                    mateStart=trChim[chimRepresent]->exons[iex][EX_G];
                    mateStrand=(uint8_t) (trChim[chimRepresent]->Str!=trChim[chimRepresent]->exons[iex][EX_iFrag]);
                };
            };

        };

        bamN+=RA->alignBAM(*trChim[itr], chimN, iTr, RA->mapGen.chrStart[trChim[itr]->Chr],  mateChr, \
                           mateStart-RA->mapGen.chrStart[(mateChr<RA->mapGen.nChrReal ? mateChr : 0)], mateStrand, alignType, \
                           NULL, P.outSAMattrOrder, RA->outBAMoneAlign+bamN, RA->outBAMoneAlignNbytes+bamN);
        bamBytesTotal+=RA->outBAMoneAlignNbytes[0]+RA->outBAMoneAlignNbytes[1];//outBAMoneAlignNbytes[1] = 0 if SE is recorded
    };

    //write all bam lines
    for (int ii=0; ii<bamN; ii++) {//output all pieces
        int tagI=-1;
        if (ii==bamIrepr) {
            tagI=bamIsuppl;
        } else if (ii==bamIsuppl) {
            tagI=bamIrepr;
        };
        if (tagI>=0) {
            bam1_t *b;
            b=bam_init1();
            bam_read1_fromArray(RA->outBAMoneAlign[tagI], b);
            uint8_t* auxp=bam_aux_get(b,"NM");
            uint32_t auxv=bam_aux2i(auxp);
            string tagSA1="SAZ"+RA->mapGen.chrName[b->core.tid]+','+to_string((uint)b->core.pos+1) +',' + ( (b->core.flag&0x10)==0 ? '+':'-') + \
                    ',' + bam_cigarString(b) + ',' + to_string((uint)b->core.qual) + ',' + to_string((uint)auxv) + ';' ;

            memcpy( (void*) (RA->outBAMoneAlign[ii]+RA->outBAMoneAlignNbytes[ii]), tagSA1.c_str(), tagSA1.size()+1);//copy string including \0 at the end
            RA->outBAMoneAlignNbytes[ii]+=tagSA1.size()+1;
             * ( (uint32*) RA->outBAMoneAlign[ii] ) = RA->outBAMoneAlignNbytes[ii]-sizeof(uint32);
	    free(b); // don't use bam_destroy1(), because bam_read1_fromArray does not allocate memory for b->data
        };

        if (P.outBAMunsorted) RA->outBAMunsorted->unsortedOneAlign(RA->outBAMoneAlign[ii], RA->outBAMoneAlignNbytes[ii], ii>0 ? 0 : bamBytesTotal);
        if (P.outBAMcoord)    RA->outBAMcoord->coordOneAlign(RA->outBAMoneAlign[ii], RA->outBAMoneAlignNbytes[ii], (RA->iReadAll<<32) );
    };

};
