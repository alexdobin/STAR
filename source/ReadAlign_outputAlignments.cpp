#include "ReadAlign.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

void ReadAlign::outputAlignments() {
  
    outBAMbytes=0;
        
    readAnnot.reset();
    
    if (mapGen.pGe.gType==101) {//temporary
        ReadAlign::spliceGraphWriteSAM();
        return;
    };    

    ReadAlign::outFilterBySJout();//sets outFilterBySJoutPass=false if read is held for the 2nd stage of outFilterBySJout
    
    if (outFilterBySJoutPass) {//otherwise align is held for the 2nd stage of outFilterBySJout
        ////////////////////////////////////
        if (unmapType<0) {//passed mappedFilter. Unmapped reads can have nTr>0

            auto nTr1 = nTr;
            auto trOut1 = trMult[0];

            if (P.pGe.transform.outYes) {
                nTr1 = alignsGenOut.alN;
                trOut1 = alignsGenOut.alMult[0];
            };

            if (nTr1>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType = -2; //not sure if this used
            } else if (nTr1==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*trOut1, Lread);
            };

            if (P.pGe.transform.outSAM && (!P.twoPass.yes || P.twoPass.pass2) ) {//transform genome only on 2nd pass
                ReadAlign::recordSJ(alignsGenOut.alN, alignsGenOut.alMult, chunkOutSJ);
            } else {
                ReadAlign::recordSJ(nTr, trMult, chunkOutSJ); //this will set mateMapped
            };            
            
            ReadAlign::alignedAnnotation();
        };

        //the operations below are both for mapped and unmapped reads
        soloRead->readBar->getCBandUMI(Read0, Qual0, readLengthOriginal, readNameExtra[0], readFilesIndex, readName);

        //transcripts: need to be run after CB/UMI are obtained to output CR/UR tags
        if ( P.quant.trSAM.yes && unmapType<0) {//Aligned.toTranscriptome output, only for mapped
            if (P.pGe.transform.outQuant) {
                quantTranscriptome(chunkTr, alignsGenOut.alN,  alignsGenOut.alMult,  alignTrAll);
            } else {
                quantTranscriptome(chunkTr, nTr, trMult,  alignTrAll);
            };
        };        
        
        soloRead->record((unmapType<0 ? nTr : 0), trMult, iReadAll, readAnnot); //need to supply nTr=0 for unmapped reads

        if (P.pGe.transform.outSAM) {
            ReadAlign::writeSAM(alignsGenOut.alN, alignsGenOut.alMult, alignsGenOut.alBest);
        } else {
            ReadAlign::writeSAM(nTr, trMult, trBest); //this will set mateMapped
        };
    };    

    if (unmapType>=0) {//unmapped reads
        statsRA.unmappedAll++; //include unmapType==4, i.e. one-mate alignments of PE reads - which may have been set in writeSAM above
        ReadAlign::outReadsUnmapped(); //uses mateMapped that was set in writeSAM above
    };    
};


///////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::recordSJ(uint64 nTrO, Transcript **trO, OutSJ *cSJ)
{//junction output for mapped reads (i.e. passed BySJout filtering)
    if (!P.outSJ.yes)
        return; //no SJ output
    
    if ( P.outSJfilterReads=="All" || nTrO==1 ) {
        uint64 sjReadStartN=cSJ->N;
        for (uint64 iTr=0; iTr<nTrO; iTr++) {//write all transcripts junctions
            outputTranscriptSJ (*(trO[iTr]), nTrO, cSJ, sjReadStartN);
        };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outFilterBySJout()
{//filtering by SJout
    outFilterBySJoutPass=true;//only false if the alignment is held for outFilterBySJoutStage. True even if unmapped
    
    if (unmapType>0 || P.outFilterBySJoutStage!=1)
        return; //unmapped, or 2nd stage
   
    for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
        for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
            if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                outFilterBySJoutPass=false;
                break;
            };
        };
        if (!outFilterBySJoutPass) 
            break;
    };
    
    if (!outFilterBySJoutPass) {//this read is held for further filtering BySJout, record fastq
        unmapType=-3; //the read is not conisdered mapped
        statsRA.readN--;
        statsRA.readBases -= readLength[0]+readLength[1];

        for (uint im=0;im<P.readNends;im++) {
            chunkOutFilterBySJoutFiles[im] << readNameMates[im] <<" "<< iReadAll <<" "<< readFilter <<" "<< readFilesIndex;
            if (!readNameExtra[im].empty())
                chunkOutFilterBySJoutFiles[im]<<" "<< readNameExtra[im];
            chunkOutFilterBySJoutFiles[im] <<"\n";
            chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutFilterBySJoutFiles[im] << "+\n";
                chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
            };
        };
    };
    
    //SJ output for all reads, including those not passed bySJout filtering. This only needs to be at the 1st stage of BySJout filtering
    ReadAlign::recordSJ(nTr, trMult, chunkOutSJ1); //this will set mateMapped
         
};

////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::writeSAM(uint64 nTrOutSAM, Transcript **trOutSAM, Transcript *trBestSAM)
{
    outBAMbytes=0;
    mateMapped[0] = mateMapped[1] = false; //mateMapped = are mates present in any of the transcripts?

    if (unmapType < 0 && outFilterBySJoutPass) {//write to SAM/BAM
        
        //////////////////////////////////////////////////////////////////////////////////
        /////////////outSAMfilter
        if (P.outSAMfilter.yes) {
            if (P.outSAMfilter.KeepOnlyAddedReferences) {
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr<mapGen.genomeInsertChrIndFirst) {
                        return;//no SAM output
                    };
                };
            } else if (P.outSAMfilter.KeepAllAddedReferences) {
                uint64 nTrOutSAM1=0;
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr>=mapGen.genomeInsertChrIndFirst) {
                        trOutSAM[nTrOutSAM1]=trOutSAM[itr];
                        trOutSAM[nTrOutSAM1]->primaryFlag=false;
                        ++nTrOutSAM1;
                    };
                };
                if (nTrOutSAM1==0) {
                   return;//no SAM output
                } else {
                    trOutSAM[0]->primaryFlag=true;
                };
                nTrOutSAM = nTrOutSAM1;
            };
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////// write SAM/BAM 
        auto nTrOutWrite=min(P.outSAMmultNmax,nTrOutSAM); //number of aligns to write to SAM/BAM files            
        
        for (uint iTr=0;iTr<nTrOutWrite;iTr++) {//write transcripts
            //mateMapped1 = true if a mate is present in this transcript
            bool mateMapped1[2]={false,false};
            mateMapped1[trOutSAM[iTr]->exons[0][EX_iFrag]]=true;
            mateMapped1[trOutSAM[iTr]->exons[trOutSAM[iTr]->nExons-1][EX_iFrag]]=true;

            if (P.outSAMbool) {//SAM output
                outBAMbytes+=outputTranscriptSAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, (uint) -1, (uint) -1, 0, -1, NULL, outSAMstream);
                if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                    outBAMbytes+= outputTranscriptSAM(*(trOutSAM[iTr]), 0, 0, (uint) -1, (uint) -1, 0, 4, mateMapped1, outSAMstream);
                };
            };

            if (P.outBAMunsorted || P.outBAMcoord) {//BAM output
                alignBAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrder,outBAMoneAlign, outBAMoneAlignNbytes);

                if (P.outBAMunsorted) {//unsorted
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite);
                    };
                    if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                        alignBAM(*trOutSAM[iTr], 0, 0, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, 4, mateMapped1, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                        for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                            outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite);
                        };
                    };
                };

                if (P.outBAMcoord) {//coordinate sorted
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (iReadAll<<32) | (iTr<<8) | trOutSAM[iTr]->exons[0][EX_iFrag] );
                    };
                };
            };
        };

        /////////////////////////////////////////////////////////////////////////////////////////////
        //////// write unmapped ends
        //TODO it's better to check all transcripts in the loop above for presence of both mates
        mateMapped[trBestSAM->exons[0][EX_iFrag]] = true;
        mateMapped[trBestSAM->exons[trBestSAM->nExons-1][EX_iFrag]] = true;

        if (P.readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) {//not readNends: this is alignment
            unmapType=4;
        };

        if (unmapType==4 && P.outSAMunmapped.within) {//output unmapped ends for single-end alignments of PE reads
            if (P.outSAMbool && !P.outSAMunmapped.keepPairs ) {
                outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
            };

            if ( P.outBAMcoord || (P.outBAMunsorted && !P.outSAMunmapped.keepPairs) ) {//BAM output
                alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                for (uint imate=0; imate<P.readNmates; imate++) {//alignBAM output is empty for mapped mate, but still need to scan through it //not readNends: this is alignment
                    if (P.outBAMunsorted && !P.outSAMunmapped.keepPairs) {
                        outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                    };
                    if (P.outBAMcoord) {//KeepPairs option does not affect for sorted BAM since we do not want multiple entries for the same unmapped read
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32);
                    };
                };
            };
        };  
        
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////write completely unmapped
    } else if (unmapType>=0 && P.outSAMunmapped.within) {//output unmapped within && unmapped read && both mates unmapped
        if (P.outBAMcoord || P.outBAMunsorted || P.quant.trSAM.bamYes) {//BAM output
            alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
            for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                if (P.outBAMunsorted) {
                    outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                };
                if (P.quant.trSAM.bamYes) {
                    outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                };
                if (P.outBAMcoord) {
                    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32);
                };
            };
        };

        if (P.outSAMbool) {//output SAM
            outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
        };
    };       
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outReadsUnmapped()
{
    if (P.outReadsUnmapped=="Fastx" ) {//output to fasta/q files. Include unmapType==4, i.e. one-mate alignments of PE reads
       for (uint im=0;im<P.readNends;im++) {
           chunkOutUnmappedReadsStream[im] << readNameMates[im]  <<" "<<im<<":"<< readFilter <<": "<< readNameExtra[im];
           if (P.readNmates>1) //not readNends: this is alignment
               chunkOutUnmappedReadsStream[im] <<" "<< int(mateMapped[0]) <<  int(mateMapped[1]);
           chunkOutUnmappedReadsStream[im] <<"\n";
           chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutUnmappedReadsStream[im] << "+\n";
                chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
            };
       };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::spliceGraphWriteSAM()
{//temporary: SAM output for SpliceGraph
    outBAMbytes=0;
    uint64 nTrOutSAM = nTr;
    if (mapGen.genomeOut.convYes) {//convert to new genome
        nTrOutSAM=0;
        for (uint32 iTr=0; iTr<nTrOutSAM; iTr++) {//convert output transcripts into new genome
            *alignsGenOut.alMult[nTrOutSAM] = *trMult[iTr];//copy information before conversion
            if (trMult[iTr]->convertGenomeCigar(*mapGen.genomeOut.g, *alignsGenOut.alMult[nTrOutSAM])) {
                ++nTrOutSAM;
                trMult[nTrOutSAM-1] = alignsGenOut.alMult[nTrOutSAM-1]; //point to new transcsript
            };
        };
        nTrOutSAM=nTrOutSAM;
    };

    for (uint iTr=0; iTr<nTrOutSAM; iTr++) {//write all transcripts            
        outBAMbytes += outputSpliceGraphSAM(*(trMult[iTr]), nTrOutSAM, iTr, outSAMstream);
    };
};

void ReadAlign::alignedAnnotation()
{
    //TODO maybe initialize readAnnot to all empty?
    //genes
    if ( P.quant.geCount.yes ) {
        if (P.pGe.transform.outQuant) {
            chunkTr->geneCountsAddAlign(alignsGenOut.alN, alignsGenOut.alMult, readAnnot.geneExonOverlap);
        } else {
            chunkTr->geneCountsAddAlign(nTr, trMult, readAnnot.geneExonOverlap);
        };        
    };
    //solo-GeneFull
    if ( P.quant.geneFull.yes ) {
        chunkTr->geneFullAlignOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull]);
    };   
    //solo-Gene
    if ( P.quant.gene.yes ) {
        chunkTr->classifyAlign(trMult, nTr, readAnnot);
    };
    //solo-GeneFull_ExonOverIntron
    if ( P.quant.geneFull_ExonOverIntron.yes ) {
        chunkTr->geneFullAlignOverlap_ExonOverIntron(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_ExonOverIntron], readAnnot.annotFeatures[SoloFeatureTypes::Gene]);
    };
    //solo-GeneFull_Ex50pAS
    if ( P.quant.geneFull_Ex50pAS.yes ) {
        chunkTr->alignExonOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_Ex50pAS]);
    };    
};