#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "SoloCommon.h"
#include <unordered_map>
#include <bitset>

#include "systemFunctions.h"


void SoloFeature::countVelocyto()
{//velocyto counting gets info from Gene counting
    time_t rawTime;

    nReadPerCB.resize(nCB);

    vector<unordered_map<uintUMI,vector<trTypeStruct>>> cuTrTypes (nCB);
    for (uint32 ii=0; ii<nCB; ii++)
        cuTrTypes[ii].reserve(readFeatSum->cbReadCount[ii] > 100 ? readFeatSum->cbReadCount[ii] : readFeatSum->cbReadCount[ii]/5); //heuristics...
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: allocated arrays" <<endl;
    
    //////////// input records
    for (int iThread=0; iThread<P.runThreadN; iThread++) {//TODO: this can be parallelized
        fstream *streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0,ios::beg);
        
        uint64 iread;
        while (*streamReads >> iread) {//until the end of file
            uintCB cb=soloFeatAll[pSolo.featureInd[SoloFeatureTypes::Gene]]->readInfo[iread].cb;
            uintUMI umi=soloFeatAll[pSolo.featureInd[SoloFeatureTypes::Gene]]->readInfo[iread].umi;
            if (cb == (uintCB)-1 || umi == (uintUMI)-1 ) {//CB and/or UMI undefined.   TODO: put a filter on CBs here, e.g. UMI threshold
                streamReads->ignore((uint32)-1, '\n');
                continue;
            };

            uint32 iCB=indCBwl[cb];
            nReadPerCB[iCB]++;//simple estimate
            
            if (cuTrTypes[iCB].count(umi)>0 && cuTrTypes[iCB][umi].empty()) {//intersection is empty, no need to load this transcript
                streamReads->ignore((uint32)-1, '\n');
                continue;
            };

            uint32 nTr;
            *streamReads >> nTr;
            vector<trTypeStruct> trT(nTr);
            for (auto & tt: trT) {
                uint32 ty;
                *streamReads >> tt.tr >> ty;
                tt.type=(uint8) ty;
            };

            if (cuTrTypes[iCB].count(umi)==0) {//1st entry for this umi
                cuTrTypes[iCB][umi]=trT;
                continue;
            };
            
            uint32 inew=0;
            vector<trTypeStruct> trT1;
            trT1.reserve(cuTrTypes[iCB][umi].size());
            
            for (uint32 iold=0; iold<cuTrTypes[iCB][umi].size(); iold++) {//intersection of old with new
                while (inew < trT.size() && cuTrTypes[iCB][umi][iold].tr>trT[inew].tr) //move through the sorted lists
                    ++inew;
                
                if (inew == trT.size() ) //end of trT reached
                    break;
                
                if (cuTrTypes[iCB][umi][iold].tr == trT[inew].tr) {//
                    trT1.push_back({trT[inew].tr, (uint8)(cuTrTypes[iCB][umi][iold].type | trT[inew].type)});
                };
            };
            cuTrTypes[iCB][umi]=trT1;//replace with intersection
        };        
    };   
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: finished input" <<endl;

    
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////// counts for  each CB
    nUMIperCB.resize(nCB,0);
    nGenePerCB.resize(nCB,0);
       
    countMatStride=4; //spliced/unspliced/ambiguous
    countCellGeneUMI.resize(nReadsMapped*countMatStride/5+16); //5 is heuristic, will be resized if needed
    countCellGeneUMIindex.resize(nCB+1);
    countCellGeneUMIindex[0]=0;
    
    for (uint32 iCB=0; iCB<nCB; iCB++) {//main collapse cycle
        map<uint32,array<uint32,3>> geneC;
        for (auto &umi : cuTrTypes[iCB]) {//cycle over UMIs
            if (umi.second.empty()) //no transcripts in the intesect
                continue;
            uint32 geneI=Trans.trGene[umi.second[0].tr];
            
            //for each transcript, are all UMIs:
            bool exonModel=false;   //purely exonic
            bool intronModel=false; //purely intronic
            bool spanModel=true; //spanning model for all UMIs
            bool mixedModel=false;  //both intronic and exonic, but not spanning
            
            for (auto &tt : umi.second) {//cycle over transcripts in this UMI
                if (Trans.trGene[tt.tr] != geneI) {//multigene
                    geneI=(uint32)-1;
                    break;
                };
                
                bitset<velocytoTypeGeneBits> gV (tt.type);
                
                mixedModel  |= ((gV.test(AlignVsTranscript::Intron) && gV.test(AlignVsTranscript::Concordant)) || gV.test(AlignVsTranscript::ExonIntron)) && !gV.test(AlignVsTranscript::ExonIntronSpan);//has exon and intron, but no span
                spanModel   &= gV.test(AlignVsTranscript::ExonIntronSpan);
                exonModel   |= gV.test(AlignVsTranscript::Concordant) && !gV.test(AlignVsTranscript::Intron) && !gV.test(AlignVsTranscript::ExonIntron);//has only-exons model                
                intronModel |= gV.test(AlignVsTranscript::Intron) && !gV.test(AlignVsTranscript::ExonIntron) && !gV.test(AlignVsTranscript::Concordant);//has only-introns model, this includes span
            };
            
            if (geneI+1==0) //multigene
                continue;
            
            if (exonModel && !intronModel && !mixedModel) {// all models are only-exons
                geneC[geneI][0]++; //spliced 
            } else if ( spanModel ||
                        ( (intronModel || mixedModel) && !exonModel)
                      ) {//all models are only-spans. Not sure if this ever happens, should be simplified. Intron or mixed, but not exon
                geneC[geneI][1]++;//unspliced
            } else {//all other combinations are mixed models, i.e. both only-exons and only-introns, only-exons and mixed
                geneC[geneI][2]++;//ambiguous
            };
            
            nUMIperCB[iCB]++;
        };
        
        countCellGeneUMIindex[iCB+1] = countCellGeneUMIindex[iCB];        
        
        if (nUMIperCB[iCB]==0) //no UMIs counted for this CB
            continue;
        
        nGenePerCB[iCB]+=geneC.size();
        readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[iCB];
        ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];
        
        if (countCellGeneUMI.size() < countCellGeneUMIindex[iCB+1] + geneC.size()*countMatStride) //allocated vector too small
            countCellGeneUMI.resize(countCellGeneUMI.size()*2);        
        
        for (auto & gg : geneC) {
            countCellGeneUMI[countCellGeneUMIindex[iCB+1]+0]=gg.first;
            for (uint32 ii=0;ii<3;ii++)
                countCellGeneUMI[countCellGeneUMIindex[iCB+1]+1+ii]=gg.second[ii];
            countCellGeneUMIindex[iCB+1] += countMatStride;
        };
    };
    
    //for now only unique mappers are considered for velocyto
    nReadPerCBtotal = nReadPerCB;
    nReadPerCBunique = nReadPerCB;
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: finished collapsing UMIs" <<endl;
    P.inOut->logMain << "RAM for solo feature "<< SoloFeatureTypes::Names[featureType] <<"\n"
                     <<  linuxProcMemory() << flush;   
};
