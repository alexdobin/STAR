#include "ReadAlign.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

void ReadAlign::outputAlignments() {

    outBAMbytes=0;

    bool mateMapped[2]={false,false};

    if (unmapType==-1) {//output transcripts

        outFilterPassed=true;

        if (P.outFilterBySJoutStage==1) {//filtering by SJout
            for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
                for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
                    if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                        outFilterPassed=false;
                        break;
                    };
                };
                if (!outFilterPassed) break;
            };
            if (!outFilterPassed) {//this read is held for further filtering BySJout, record fastq
                unmapType=-3; //the read is not conisddred unmapped
                statsRA.readN--;
                statsRA.readBases -= readLength[0]+readLength[1];

//                 if (P.runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutFilterBySJout);
                for (uint im=0;im<P.readNmates;im++) {
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
//                 if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutFilterBySJout);
            };
        };

        if (P.outSJfilterReads=="All" || nTr==1) {
            uint sjReadStartN=chunkOutSJ1->N;
            for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ1, sjReadStartN);
            };
        };

        if (outFilterPassed) {
            uint nTrOut=nTr; //number of aligns to output
            bool outSAMfilterYes=true;
            if (P.outSAMfilter.yes)
            {
                if (P.outSAMfilter.KeepOnlyAddedReferences)
                {
                    for (uint itr=0;itr<nTr;itr++)
                    {//check if transcripts map to chr other than added references
                        if (trMult[itr]->Chr<mapGen.genomeInsertChrIndFirst)
                        {
                            outSAMfilterYes=false;
                            break;
                        };
                    };
                } else if (P.outSAMfilter.KeepAllAddedReferences)
                {
                    nTrOut=0;
                    for (uint itr=0;itr<nTr;itr++)
                    {//check if transcripts map to chr other than added references
                        if (trMult[itr]->Chr>=mapGen.genomeInsertChrIndFirst)
                        {
                            trMult[nTrOut]=trMult[itr];
                            trMult[nTrOut]->primaryFlag=false;
                            ++nTrOut;
                        };
                    };
                    if (nTrOut==0)
                    {
                        outSAMfilterYes=false;
                    } else
                    {
                        trMult[0]->primaryFlag=true;
                    };
                };
            };
            if (nTr>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType=-1;
            } else if (nTr==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*(trMult[0]),Lread);
                unmapType=-2;
            } else {//cannot be
                ostringstream errOut;
                errOut  << "EXITING because of a BUG: nTr=0 in outputAlignments.cpp";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };

            nTrOut=min(P.outSAMmultNmax,nTrOut); //number of to write to SAM/BAM files

            for (uint iTr=0;iTr<nTrOut;iTr++)
            {//write all transcripts
                //mate mapped = true if a mate was present in one of the trancsripts
                //mateMapped[trMult[iTr]->exons[0][EX_iFrag]]=true;
                //mateMapped[trMult[iTr]->exons[trMult[iTr]->nExons-1][EX_iFrag]]=true;

                //mateMapped1 = true if a mate is present in this transcript
                bool mateMapped1[2]={false,false};
                mateMapped1[trMult[iTr]->exons[0][EX_iFrag]]=true;
                mateMapped1[trMult[iTr]->exons[trMult[iTr]->nExons-1][EX_iFrag]]=true;

                if (P.outSAMbool && outSAMfilterYes)
                {//SAM output
                    outBAMbytes+=outputTranscriptSAM(*(trMult[iTr]), nTr, iTr, (uint) -1, (uint) -1, 0, -1, NULL, outSAMstream);
                    if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) )
                    {//keep pairs && paired reads && one of the mates not mapped in this transcript
                        outBAMbytes+= outputTranscriptSAM(*(trMult[iTr]), 0, 0, (uint) -1, (uint) -1, 0, 4, mateMapped1, outSAMstream);
                    };
                };

                if ((P.outBAMunsorted || P.outBAMcoord) && outSAMfilterYes)
                {//BAM output
                    alignBAM(*(trMult[iTr]), nTr, iTr, mapGen.chrStart[trMult[iTr]->Chr], (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrder,outBAMoneAlign, outBAMoneAlignNbytes);

                    if (P.outBAMunsorted)
                    {//unsorted
                        for (uint imate=0; imate<P.readNmates; imate++)
                        {//output each mate
                            outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOut);
                        };
                        if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) )
                        {//keep pairs && paired reads && one of the mates not mapped in this transcript
                            alignBAM(*trMult[iTr], 0, 0, mapGen.chrStart[trMult[iTr]->Chr], (uint) -1, (uint) -1, 0, 4, mateMapped1, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                            for (uint imate=0; imate<P.readNmates; imate++)
                            {//output each mate
                                outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOut);
                            };
                        };
                    };

                    if (P.outBAMcoord)
                    {//coordinate sorted
                        for (uint imate=0; imate<P.readNmates; imate++)
                        {//output each mate
                            outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (iReadAll<<32) | (iTr<<8) | trMult[iTr]->exons[0][EX_iFrag] );
                        };
                    };

                };
            };

            mateMapped[trBest->exons[0][EX_iFrag]]=true;
            mateMapped[trBest->exons[trBest->nExons-1][EX_iFrag]]=true;            
            
            if (P.readNmates>1 && !(mateMapped[0] && mateMapped[1]) )
            {
                unmapType=4;
            };

            if (unmapType==4 && P.outSAMunmapped.yes)
            {//output unmapped end for single-end alignments
                if (P.outSAMbool && !P.outSAMunmapped.keepPairs && outSAMfilterYes)
                {
                    outBAMbytes+= outputTranscriptSAM(*trBest, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
                };

                if ( (P.outBAMcoord || (P.outBAMunsorted && !P.outSAMunmapped.keepPairs) ) && outSAMfilterYes)
                {//BAM output
                    alignBAM(*trBest, 0, 0, mapGen.chrStart[trBest->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                    for (uint imate=0; imate<P.readNmates; imate++)
                    {//alignBAM output is empty for mapped mate, but still need to scan through it
                        if (P.outBAMunsorted && !P.outSAMunmapped.keepPairs)
                        {
                            outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                        };
                        if (P.outBAMcoord)
                        {//KeepPairs option does not affect for sorted BAM since we do not want multiple entries for the same unmapped read
                            outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll);
                        };
                    };
                };
            };

            if (P.outSJfilterReads=="All" || nTr==1)
            {
                uint sjReadStartN=chunkOutSJ->N;
                for (uint iTr=0;iTr<nTr;iTr++)
                {//write all transcripts
                    outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ, sjReadStartN);
                };
            };

            if ( P.quant.geCount.yes )
            {
                chunkTr->geneCountsAddAlign(nTr, trMult);
            };

            if ( P.quant.trSAM.yes )
            {//NOTE: the transcripts are changed by this function (soft-clipping extended), cannot be reused
                quantTranscriptome(chunkTr, nTrOut, trMult,  alignTrAll);
            };
        };
    };

    if (unmapType>=0)
    {
        statsRA.unmappedAll++;
    };

    if ( P.outSAMunmapped.within && unmapType>=0 && unmapType<4 ) {//output unmapped within && unmapped read && both mates unmapped
        if (P.outBAMcoord || P.outBAMunsorted || P.quant.trSAM.yes)
        {//BAM output
            alignBAM(*trBest, 0, 0, mapGen.chrStart[trBest->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
            for (uint imate=0; imate<P.readNmates; imate++)
            {//output each mate
                if (P.outBAMunsorted)
                {
                    outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                };
                if (P.quant.trSAM.yes)
                {
                    outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
                };
                if (P.outBAMcoord)
                {
                    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll);
                };
            };
        };

        if (P.outSAMbool)
        {//output SAM
            outBAMbytes+= outputTranscriptSAM(*trBest, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
        };
    };
    if (unmapType>=0 && P.outReadsUnmapped=="Fastx" ){//output to fasta/q files
           for (uint im=0;im<P.readNmates;im++) {
               chunkOutUnmappedReadsStream[im] << readNameMates[im];
               if (P.readNmates>1) chunkOutUnmappedReadsStream[im] <<"\t"<< int(mateMapped[0]) <<  int(mateMapped[1]);
               chunkOutUnmappedReadsStream[im] <<"\n";
               chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
                if (readFileType==2) {//fastq
                    chunkOutUnmappedReadsStream[im] << "+\n";
                    chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
                };
           };
    };
};



