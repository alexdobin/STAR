#include "ReadAlign.h"
#include "BAMfunctions.h"

void ReadAlign::chimericDetectionOldOutput() {
    
    if (!chimRecord) {
        return;
    };
    
    chimN=2; //this  is hard-coded for now
    //re-calculate the score for chimeric transcripts
    trChim[0].alignScore(Read1, mapGen.G, P);
    trChim[1].alignScore(Read1, mapGen.G, P);

    int chimRepresent=-999, chimType=0;
    if (trChim[0].exons[0][EX_iFrag]!=trChim[0].exons[trChim[0].nExons-1][EX_iFrag]) {//tr0 has both mates
        chimRepresent = 0;
        chimType = 1;
        trChim[0].primaryFlag=true;//paired portion is primary
        trChim[1].primaryFlag=false;
    } else if (trChim[1].exons[0][EX_iFrag]!=trChim[1].exons[trChim[1].nExons-1][EX_iFrag]) {//tr1 has both mates
        chimRepresent = 1;
        chimType = 1;
        trChim[1].primaryFlag=true;//paired portion is primary
        trChim[0].primaryFlag=false;
    } else if (trChim[0].exons[0][EX_iFrag]!=trChim[1].exons[0][EX_iFrag]) {//tr0 and tr1 are single different mates
        chimRepresent = -1;
        chimType = 2;
        trChim[0].primaryFlag=true;
        trChim[1].primaryFlag=true;
    } else  {//two chimeric segments are on the same mate - this can only happen for single-end reads
        chimRepresent = (trChim[0].maxScore > trChim[1].maxScore) ? 0 : 1;
        chimType = 3;
        trChim[chimRepresent].primaryFlag=true;
        trChim[1-chimRepresent].primaryFlag=false;
    };

    if (P.pCh.out.bam) {//BAM output
        int alignType, bamN=0, bamIsuppl=-1, bamIrepr=-1;
        uint bamBytesTotal=0;//estimate of the total size of all bam records, for output buffering
        uint mateChr,mateStart;
        uint8_t mateStrand;
        for (uint itr=0;itr<chimN;itr++) {//generate bam for all chimeric pieces
            if (chimType==2) {//PE, encompassing
                mateChr=trChim[1-itr].Chr;
                mateStart=trChim[1-itr].exons[0][EX_G];
                mateStrand=(uint8_t) (trChim[1-itr].Str!=trChim[1-itr].exons[0][EX_iFrag]);
                alignType=-10;
            } else {//spanning chimeric alignment, could be PE or SE
                mateChr=-1;mateStart=-1;mateStrand=0;//no need fot mate info unless this is the supplementary alignment
                if (chimRepresent==(int)itr) {
                    alignType=-10; //this is representative part of chimeric alignment, record is as normal; if encompassing chimeric junction, both are recorded as normal
                    bamIrepr=bamN;
                    if (trChim[itr].exons[0][EX_iFrag]!=trChim[1-itr].exons[0][EX_iFrag]) {//the next mate is chimerically split
                        ++bamIrepr;
                    };
//                     if (chimType==3) {
//                         bamIrepr=bamN;
//                     } else if (trChim[itr].exons[0][EX_iFrag]==trChim[1-itr].exons[0][EX_iFrag]) {
//                         
//                    };
//                     bamIrepr=( (itr%2)==(trChim[itr].Str) && chimType!=3) ? bamN+1 : bamN;//this is the mate that is chimerically split
                } else {//"supplementary" chimeric segment
                    alignType=P.pCh.out.bamHardClip ? ( ( itr%2==trChim[itr].Str ) ? -12 : -11) : -13 ; //right:left chimeric junction
                    bamIsuppl=bamN;
                    if (chimType==1) {//PE alignment, need mate info for the suppl
                        uint iex=0;
                        for (;iex<trChim[chimRepresent].nExons-1;iex++) {
                            if (trChim[chimRepresent].exons[iex][EX_iFrag]!=trChim[itr].exons[0][EX_iFrag]) {
                                break;
                            };
                        };
                        mateChr=trChim[chimRepresent].Chr;
                        mateStart=trChim[chimRepresent].exons[iex][EX_G];
                        mateStrand=(uint8_t) (trChim[chimRepresent].Str!=trChim[chimRepresent].exons[iex][EX_iFrag]);
                    };
                };

            };

            bamN+=alignBAM(trChim[itr], 1, 0, mapGen.chrStart[trChim[itr].Chr],  mateChr, mateStart-mapGen.chrStart[(mateChr<mapGen.nChrReal ? mateChr : 0)], mateStrand, \
                            alignType, NULL, P.outSAMattrOrder, outBAMoneAlign+bamN, outBAMoneAlignNbytes+bamN);
            bamBytesTotal+=outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1];//outBAMoneAlignNbytes[1] = 0 if SE is recorded
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
                bam_read1_fromArray(outBAMoneAlign[tagI], b);
                uint8_t* auxp=bam_aux_get(b,"NM");
                uint32_t auxv=bam_aux2i(auxp);
                string tagSA1="SAZ"+mapGen.chrName[b->core.tid]+','+to_string((uint)b->core.pos+1) +',' + ( (b->core.flag&0x10)==0 ? '+':'-') + \
                        ',' + bam_cigarString(b) + ',' + to_string((uint)b->core.qual) + ',' + to_string((uint)auxv) + ';' ;

                 memcpy( (void*) (outBAMoneAlign[ii]+outBAMoneAlignNbytes[ii]), tagSA1.c_str(), tagSA1.size()+1);//copy string including \0 at the end
                 outBAMoneAlignNbytes[ii]+=tagSA1.size()+1;
                 * ( (uint32*) outBAMoneAlign[ii] ) = outBAMoneAlignNbytes[ii]-sizeof(uint32);
            };

            if (P.outBAMunsorted) outBAMunsorted->unsortedOneAlign(outBAMoneAlign[ii], outBAMoneAlignNbytes[ii], ii>0 ? 0 : bamBytesTotal);
            if (P.outBAMcoord)    outBAMcoord->coordOneAlign(outBAMoneAlign[ii], outBAMoneAlignNbytes[ii], (iReadAll<<32) );
        };
    };

    if (P.pCh.out.samOld) {
        for (uint iTr=0;iTr<chimN;iTr++) 
        {//write all chimeric pieces to Chimeric.out.sam/junction
            if (P.readNmates==2) {//PE: need mate info
                uint iex=0;
                if ( trChim[1-iTr].exons[0][EX_iFrag] != trChim[1-iTr].exons[trChim[1-iTr].nExons-1][EX_iFrag] )
                {//the other segment has 2 mates, need to find the opposite mate
                    for (;iex<trChim[1-iTr].nExons;iex++) {
                        if (trChim[1-iTr].exons[iex][EX_iFrag]!=trChim[iTr].exons[0][EX_iFrag]) {
                            break;
                        };
                    };
                };

                uint mateChr=trChim[1-iTr].Chr;
                uint mateStart=trChim[1-iTr].exons[iex][EX_G];
                char mateStrand=(char) (trChim[1-iTr].Str!=trChim[1-iTr].exons[iex][EX_iFrag]);

                outputTranscriptSAM(trChim[iTr], chimN, iTr, mateChr, mateStart, mateStrand, -1, NULL, &chunkOutChimSAM);
            } else 
            {
                outputTranscriptSAM(trChim[iTr], chimN, iTr, -1, -1, -1, -1, NULL, &chunkOutChimSAM);
            };
        };
    };
    
    if (P.pCh.out.junctions) {
        //junction + SAMp
        *chunkOutChimJunction << mapGen.chrName[trChim[0].Chr] <<"\t"<< chimJ0 - mapGen.chrStart[trChim[0].Chr]+1 <<"\t"<< (trChim[0].Str==0 ? "+":"-") \
                <<"\t"<< mapGen.chrName[trChim[1].Chr] <<"\t"<< chimJ1 - mapGen.chrStart[trChim[1].Chr]+1 <<"\t"<< (trChim[1].Str==0 ? "+":"-") \
                <<"\t"<< chimMotif <<"\t"<< chimRepeat0  <<"\t"<< chimRepeat1 <<"\t"<< readName+1 \
                <<"\t"<< trChim[0].exons[0][EX_G] - mapGen.chrStart[trChim[0].Chr]+1 <<"\t"<< outputTranscriptCIGARp(trChim[0]) \
                <<"\t"<< trChim[1].exons[0][EX_G] - mapGen.chrStart[trChim[1].Chr]+1 <<"\t"<<  outputTranscriptCIGARp(trChim[1]);
        if (P.outSAMattrPresent.RG)
            *chunkOutChimJunction <<"\t"<< P.outSAMattrRG.at(readFilesIndex);
        *chunkOutChimJunction <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
    };
    
    return;
};