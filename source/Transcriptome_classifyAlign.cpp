#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>

int alignToTranscript(Transcript &aG, uint trS1, uint32 *exSE1, uint16 exN1) 
{    
    bool alignIntronic      =false;
    bool alignExonic        =false;
    bool alignSpansExonIntr =false;
    bool alignSJconcordant  =true;
    
    //we assume that align is fully contained in the transcript, i.e. alignStart>=trStart, alignEnd<=trEnd
    //find exon that overlaps beginning of the read
    
    //TODO
    //iab=0;
    //distTSS=aG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1];
    //iab=aG.nExons-1;
    //distTTS=trLen[tr1]-(gthaG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1]+aG.exons[iab][EX_L]);
    //if trStr1==2: TTS=trLen[tr1]-TTS, TSS=trLen[tr1]-TSS
    
    for (uint32 iab=0, ex1=0, bS=0, bE=0, eE=0, enS=0;
                iab<aG.nExons; iab++) {//scan through all blocks of the align

        uint64 bEprev=bE;

        bS=(uint32) (aG.exons[iab][EX_G]-trS1);//block start
        bE=bS+aG.exons[iab][EX_L]-1;//block end

        if (iab==0 || aG.canonSJ[iab-1]==-3) {//start of alig, or jump to another mate
            ex1=binarySearch1<uint32>(bS, exSE1, 2*exN1) / 2;// alignStart>=ex1start
        } else if (aG.canonSJ[iab-1]>=0) {//splice junction
            if (bEprev == eE && bS == enS) {//eE and enS are still from the old ex1
                ++ex1; //junction agrees
            } else {
                alignSJconcordant = false;
                ex1=binarySearch1<uint32>(bS, exSE1, 2*exN1) / 2;// alignStart>=ex1start
            };
        };

        //uint64 eS=exSE1[2*ex1];
        eE  = exSE1[2*ex1+1];
        enS = ex1+1<exN1 ? exSE1[2*(ex1+1)] : 0;//next exon start

        if (bS <= eE) {//block starts in the ex1 exon
            if (bE > eE) {
                alignSpansExonIntr = true;
                //break;//if ex/in span is detected, no need to check anything else - no true : might still have non-concordant junction
            };
            alignExonic = true;
        } else {//block starts in the intron
            if (bE >= enS) {//block overlaps next exon
                alignSpansExonIntr = true;
                //break;//if ex/in span is detected, no need to check anything else
            };
            alignIntronic = true;
        };
    };//cycle over align blocks

    if (!alignSJconcordant) //if align has a junction, it's always checked for concordance
        return -1;          //even for exon/intron aligns, sjs have to be concordant, otherwise align is not consistent with this transcript model

    if (alignSpansExonIntr) {
        return AlignVsTranscript::ExonIntronSpan; //align spans exon/intron boundary
    } else if (!alignIntronic) {//align is purely exonic
        return AlignVsTranscript::Concordant; //align is fully exonic and concordant
    } else {//align has introns
        if (alignExonic) {
            return AlignVsTranscript::ExonIntron; //mixed exonic/intronic align, but no span
        } else {
            return AlignVsTranscript::Intron; //purely intronic align
        };
    };
};

int alignToTranscriptMinOverlap(Transcript &aG, uint trS1, uint32 *exSE1, uint16 exN1, uint32 minOverlapMinusOne) 
{
    bool alignIntronic      =false;
    bool alignExonic        =false;
    bool alignSpansExonIntr =false;
    bool alignSJconcordant  =true;
    
    //we assume that align is fully contained in the transcript, i.e. alignStart>=trStart, alignEnd<=trEnd
    //find exon that overlaps beginning of the read
    
    for (uint32 iab=0, ex1=0, bS=0, bE=0, eS=0, eE=0, enS=0; 
                iab<aG.nExons; iab++) {//scan through all blocks of the align
              
        bS=(uint32) (aG.exons[iab][EX_G]-trS1);//block start
        bE=bSt+aG.exons[iab][EX_L]-1;//block end
            
        if (iab==0 || aG.canonSJ[iab-1]==-3 || aG.canonSJ[iab-1]>=0) {//start of align, or jump to another mate, or junction
            ex1=binarySearch1<uint32>(bS, exSE1, 2*exN1) / 2;// alignStart>=ex1start            
        };

        eS  = exSE1[2*ex1];
        eE  = exSE1[2*ex1+1];
        enS = ex1+1<exN1 ? exSE1[2*(ex1+1)] : 0;//next exon start

        uint32 bStype=0;
        if (bS <= eE-minOverlapMinusOne && bS>=eS+minOverlapMinusOne) {
            bStype=1;//bS exonic
        } else if (bS>eE+minOverlapMinusOne && bS<enS-minOverlapMinusOne) {
            bStype=2;//bS intronic
        };//otherwise it's undefined, to close to exon end
        
        uint32 bEtype=0;
        if (bE >= enS+minOverlapMinusOne) {
            bEtype=1;//bE exonic
        } else if (bE<eE-minOverlapMinusOne) {
            bEtype=2;//bE intronic
        };//otherwise it's undefined, to close to exon end        
        
        if (bStype==0)
            bStype=bEtype
        
        if (bS <= eE) {//block starts in the ex1 exon
            if (bE > eE) {
                alignSpansExonIntr = true;
            };
            alignExonic = true;
        } else {//block starts in the intron
            if (bE >= enS) {//block overlaps next exon
                alignSpansExonIntr = true;
            };
            
            if (enS-eE>1000000) {//if intron is too large, do not mark this align as intronic. This is to match velocyto.py logic. TODO: make it a parameter
                return -1;//no intronic match since the intron is too large
            };

            alignIntronic = true;
        };
        
        if (aG.sjYes && (alignIntronic || alignSpansExonIntr)) {//spliced align cannot overlap an intron
            alignSJconcordant=false;
            break;
        };
    };//cycle over align blocks
    
    if (!alignSJconcordant) //if align has a junction, it's always checked for concordance
        return -1;          //even for exon/intron aligns, sjs have to be concordant, otherwise align is not consistent with this transcript model
    
    if (alignSpansExonIntr) {
        return AlignVsTranscript::ExonIntronSpan; //align spans exon/intron boundary
    } else if (!alignIntronic) {//align is purely exonic
        return AlignVsTranscript::Concordant; //align is fully exonic and concordant
     } else {//align has introns
        if (alignExonic) {
            return AlignVsTranscript::ExonIntron; //mixed exonic/intronic align, but no span
        } else {
            return AlignVsTranscript::Intron; //purely intronic align
        };
    };
};

void Transcriptome::classifyAlign (Transcript **alignG, uint64 nAlignG, ReadAnnotations &readAnnot) 
{
    readAnnot.transcriptConcordant={};
    readAnnot.geneConcordant={};

    //array<bool,AlignVsTranscript::N> reAnn={false};
    uint32 reGe=(uint32)-2;//so that the first gene can be recorded
    std::bitset<velocytoTypeGeneBits> reAnn; //initialized to 0 (false)
       
    for (uint iag=0; iag<nAlignG; iag++) {
        
        Transcript &aG=*alignG[iag];

        //binary search through transcript starts
        uint32 tr1=binarySearch1a<uint>(aG.exons[0][EX_G], trS, nTr);//tr1 has the maximum transcript start such that it is still <= align start
        if (tr1==(uint32) -1) 
            continue; //this alignment is outside of range of all transcripts

        uint aGend=aG.exons[aG.nExons-1][EX_G]+aG.exons[aG.nExons-1][EX_L]-1;

        ++tr1;
        do {//cycle back through all the transcripts
            --tr1;
            if ( aGend>trE[tr1] ||
                 (P.pSolo.strand >= 0 && (trStr[tr1]==1 ? aG.Str : 1-aG.Str) != (uint32)P.pSolo.strand) ) //!(this transcript contains the read, and has correct strand)
                     continue;
                 
            int aStatus=alignToTranscript(aG, trS[tr1], exSE+2*trExI[tr1], trExN[tr1]);
            
            if (aStatus==AlignVsTranscript::Concordant) {//align conforms with the transcript

                //TODO!!!FIX THIS for Quant!!!
                //uint32 distTTS=trLen[tr1]-(aTall[nAtr].exons[aTall[nAtr].nExons-1][EX_G] + aTall[nAtr].exons[aTall[nAtr].nExons-1][EX_L]);
                //readAnnot.transcriptConcordant.push_back({tr1,distTTS});

                readAnnot.geneConcordant.insert(trGene[tr1]);//genes for all alignments
                aG.alignGenes.insert(trGene[tr1]);//genes for each alignment
            };
            
            if (P.pSolo.featureYes[SoloFeatureTypes::Velocyto] && nAlignG==1) {//another calculation for velocyto with minOverlapMinusOne=6
                
                //6 is the hard code minOverlapMinusOne, to agree with velocyto's MIN_FLANK=5
                aStatus=alignToTranscriptMinOverlap(aG, trS[tr1], exSE+2*trExI[tr1], trExN[tr1], 6);
                
                if (aStatus<0)
                    continue; //align is not concordant with this transcript because one of the align junctions does not match transcript junction

                //calculate reAnn
                if (reGe!=(uint32)-1) {//not multi-mapper

                    if (reGe==(uint32)-2) //first gene
                        reGe=trGene[tr1];

                    if (reGe!=trGene[tr1]) {//
                        reGe=(uint32)-1; //marks multi-gene align
                    } else {
                        if (aStatus!=AlignVsTranscript::ExonIntronSpan) {
                            reAnn.set(AlignVsTranscript::ExonIntronSpan, true);//meaning of this bit is NoExonIntronSpan
                            reAnn.set(aStatus, true);
                        };
                    };
                };
            };
        } while (trEmax[tr1]>=aGend && tr1>0);
    };
    
    //velocyto logic
    readAnnot.geneVelocyto[0]=(reGe+2==0 ? (uint32)-1 : reGe);//-2 marks no gene, convert to -1 which marks either no gene or multigene - no output     
    readAnnot.geneVelocyto[1]=reAnn.to_ulong();
    
    
//     if (reAnn[AlignVsTranscript::ExonIntronSpan]) {
//         readAnnot.geneVelocyto[1]=0;
//     } else if (reAnn[AlignVsTranscript::Concordant] || reAnn[AlignVsTranscript::Exon]) {//at least one model is exonic
//         if (!reAnn[AlignVsTranscript::Intron] && !reAnn[AlignVsTranscript::ExonIntron]) {
//             readAnnot.geneVelocyto[1]=1; //spliced
//         } else {
//             readAnnot.geneVelocyto[1]=3; //ambiguous <= exonic * ( intronic || mixed)
//         };
//     } else {//all other combinations are unspliced, as they do not contain a single purely exonic model
//         readAnnot.geneVelocyto[1]=2;
//     };
    
//     if (reGe==(uint32)-1) {//multi-gene
//         readAnnot.geneVelocyto[1]=0;//value does not matter, it will not be recorded
//     } else if (reAnn[AlignVsTranscript::ExonIntronSpan]) {
//         readAnnot.geneVelocyto[1]=2; //unspliced
//     } else if (reAnn[AlignVsTranscript::Concordant] || reAnn[AlignVsTranscript::Exon]) {//at least one model is exonic
//         if (!reAnn[AlignVsTranscript::Intron] && !reAnn[AlignVsTranscript::ExonIntron]) {
//             readAnnot.geneVelocyto[1]=1; //spliced
//         } else {
//             readAnnot.geneVelocyto[1]=3; //ambiguous <= exonic * ( intronic || mixed)
//         };
//     } else {//all other combinations are unspliced, as they do not contain a single purely exonic model
//         readAnnot.geneVelocyto[1]=2;
//     };
};
