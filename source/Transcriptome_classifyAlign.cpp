#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"

int alignToTranscript(Transcript &aG, uint trS1, uint8 trStr1, uint32 *exSE1, uint32 *exLenCum1, uint16 exN1) 
{
    //returned values: 1: align fully agrees with transcript, including splices
    //                 2: align is fully exonic, but not concordant
    //                 3: align is fully intronic
    //                 4: align has blocks mapping to exons and introns, but no spanning
    //                 5: align spans exon/intron boundary
    

    bool alignIntronic      =false;
    bool alignExonic        =false;
    bool alignSpansExonIntr =false;
    bool alignSJconcordant  =true;
    
    //we assume that align is fully contained in the transcript, i.e. alignStart>=trStart, alignEnd<=trEnd
    //find exon that overlaps beginning of the read
    //uint32 g1=aG.exons[0][EX_G]-trS1;//start of the align
    //uint32 ex1=binarySearch1<uint32>(g1, exSE1, 2*exN1) / 2;// ex1start<=alignStart
    
    //TODO
    //iab=0;
    //distTSS=aG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1];
    //iab=aG.nExons-1;
    //distTTS=trLen[tr1]-(gthaG.exons[iab][EX_G]-trS1-exSE1[2*ex1]+exLenCum1[ex1]+aG.exons[iab][EX_L]);
    //if trStr1==2: TTS=trLen[tr1]-TTS, TSS=trLen[tr1]-TSS
    
    //aG.canonSJ[aG.nExons-1]=-999; //marks the last block
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
                break;//if ex/in span is detected, no need to check anything else
            };
            alignExonic = true;
        } else {//block starts in the intron
            if (bE >= enS) {//block overlaps next exon
                alignSpansExonIntr = true;
                break;//if ex/in span is detected, no need to check anything else
            };
            alignIntronic = true;
        };
    };
    
    if (alignSpansExonIntr) {
        return AlignVsTranscript::ExonIntronSpan; //align spans exon/intron boundary
    } else if (!alignIntronic) {//align is purely exonic
        if (alignSJconcordant) {
            return AlignVsTranscript::Concordant; //align is concordant with the transcript, i.e. fully agrees with transcript, including splices
        } else {
            return AlignVsTranscript::Exon; //align is fully exonic, but not concordant
        };
    } else {//align has introns
        if (alignExonic) {
            return AlignVsTranscript::ExonIntron; //mixed exonic/intronic align
        } else {
            return AlignVsTranscript::Intron; //purely intronic align
        };
    };

    return (uint32)-1; //this should not happen
};

void Transcriptome::classifyAlign (Transcript **alignG, uint64 nAlignG, ReadAnnotations &readAnnot) 
{
    readAnnot.transcriptConcordant={};
    readAnnot.geneConcordant={};

    array<bool,AlignVsTranscript::N> reAnn={false};
    uint32 reGe=(uint32)-1;
    
    for (uint iag=0; iag<nAlignG; iag++) {
        
        Transcript &aG=*alignG[iag];

        //binary search through transcript starts
        uint32 tr1=binarySearch1a<uint>(aG.exons[0][EX_G], trS, nTr);//tr1 has the maximum transcript start such that it is still <= align start
        if (tr1==(uint32) -1) 
            continue; //this alignment is outside of range of all transcripts

        uint aGend=aG.exons[aG.nExons-1][EX_G];

        ++tr1;
        do {//cycle back through all the transcripts
            --tr1;
            if ( aGend>trE[tr1] ||
                 (P.pSolo.strand >= 0 && (trStr[tr1]==1 ? aG.Str : 1-aG.Str) != (uint32)P.pSolo.strand) ) //this transcript contains the read, and has correct strand
                     continue;
                 
            int aStatus=alignToTranscript(aG, trS[tr1], trStr[tr1], exSE+2*trExI[tr1], exLenCum+trExI[tr1], trExN[tr1]);
            if (aStatus==AlignVsTranscript::Concordant) {//align conforms with the transcript

                //TODO!!!FIX THIS//uint32 distTTS=trLen[tr1]-(aTall[nAtr].exons[aTall[nAtr].nExons-1][EX_G] + aTall[nAtr].exons[aTall[nAtr].nExons-1][EX_L]);
                //readAnnot.transcriptConcordant.push_back({tr1,distTTS});

                readAnnot.geneConcordant.insert(trGene[tr1]);//genes for all alignments
                aG.alignGenes.insert(trGene[tr1]);//genes for each alignment
            };
                       
            if (reGe==(uint32)-1)
                reGe=trGene[tr1];
            if (reGe!=trGene[tr1])
                reGe=(uint32)-1; //marks multi-gene align
            
            reAnn[aStatus]=true;
            
        } while (trEmax[tr1]>=aGend && tr1>0);
    };
    
    //velocyto logic
    readAnnot.geneVelocyto[0]=reGe;
    if (reGe==(uint32)-1) {//multi-gene
        readAnnot.geneVelocyto[1]=0;
    } else if (reAnn[AlignVsTranscript::ExonIntronSpan]) {
        readAnnot.geneVelocyto[1]=2; //unspliced
    } else if (reAnn[AlignVsTranscript::Concordant] || reAnn[AlignVsTranscript::Exon]) {//at least one model is exonic
        if (!reAnn[AlignVsTranscript::Intron] && !reAnn[AlignVsTranscript::ExonIntron]) {
            readAnnot.geneVelocyto[1]=1; //spliced
        } else {
            readAnnot.geneVelocyto[1]=3; //ambiguous <= exonic * ( intronic || mixed)
        };
    } else {//all other combinations are unspliced, as they do not contain a single purely exonic model
        readAnnot.geneVelocyto[1]=2;
    };
};
