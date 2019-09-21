//
//  GTF_SuperTranscript.cpp
//  
//
//  Created by Fahimeh Mirhaj on 5/30/19.
//

#include "GTF.h"
#include "streamFuns.h"


void GTF::superTranscript() {
    /*
     Finding intervals from given gtf class structure.
     Changing the genome (given sequence) coordinates to "CondensedGene" coordinates respectively.
     assumption: the gtf data is sorted by strands and then by coordinates.
     */
    if (P.pGe.gType!="Transcriptome" && P.pGe.gType!="SuperTranscriptome")
        return;
    
    for(uint64 i = 0; i < exonLoci.size(); i++) {//convert (-)strand exons coordinates
        uint64 transId = exonLoci[i][exT];
        if(transcriptStrand[transId] == 2) {
            uint64 temp = exonLoci[i][exS];
            exonLoci[i][exS] = 2 * genome.nGenome - 1 - exonLoci[i][exE];
            exonLoci[i][exE] = 2 * genome.nGenome - 1 - temp;
        }
    }
    
    //sort by exon start position
    sort(exonLoci.begin(), exonLoci.end(),
         [](const array<uint64, exL>& e1, const array<uint64, exL>& e2) {
             return e1[exS] < e2[exS];
         });
    
    uint64 gapValue = exonLoci[0][exS];
    pair<uint64,uint64> curr = {exonLoci[0][exS], exonLoci[0][exE]}; //intervals[0]
    exonLoci[0][exS] = 0;
    exonLoci[0][exE] -= gapValue;
    vector<pair<uint64,uint64>> mergedIntervals;
    
    for(uint64 i = 1; i < exonLoci.size(); ++i) {
        if(exonLoci[i][exS] <= curr.second+1) {
            curr.second = max(curr.second, exonLoci[i][exE]);
        } else {
            gapValue += exonLoci[i][exS] - curr.second - 1;
            mergedIntervals.push_back(curr);
            curr = {exonLoci[i][exS], exonLoci[i][exE]};
        }
        exonLoci[i][exS] -= gapValue;
        exonLoci[i][exE] -= gapValue;
    }
    mergedIntervals.push_back(curr);
    
    //Condensed genome (a.k.a SuperTranscriptome) sequence
    for(const auto& p: mergedIntervals) {
        for(uint64 j = p.first; j <= p.second; ++j) {
            sequenceOfCondensedGenome.push_back((uint8) genome.G[j]);
        };
    };
    P.inOut->logMain << "SuperTranscriptome (condensed) genome length = " << sequenceOfCondensedGenome.size() <<endl;
    
    //Sort rows of exonLoci based on transId and startInterval
    //Find super and normal transcript intervals and sequences
    transcriptIntervals.resize(transcriptID.size(), {(uint64)-1, 0}); //normal transcript
    for(uint64 i = 0 ; i < exonLoci.size(); i++) {//start-end of the normal transcripts in the condensed genome coordinates
        uint64 transId = exonLoci[i][exT];
        transcriptIntervals[transId].first = min(transcriptIntervals[transId].first, exonLoci[i][exS]);
        transcriptIntervals[transId].second = max(transcriptIntervals[transId].second, exonLoci[i][exE]);
    };
    
    vector<pair<uint64, uint64>> tempTranscriptIntervals = transcriptIntervals;//sorted transript intervals
    sort(tempTranscriptIntervals.begin(), tempTranscriptIntervals.end(),
         [](const pair<uint64, uint64>& t1, const pair<uint64, uint64>& t2) {
             return t1.first < t2.first;
         });
    curr = tempTranscriptIntervals.front();
    for(const auto& p: tempTranscriptIntervals) {//if transcript intervals overlap by >=1 base, merge them into SuperTranscripts
        if(p.first <= curr.second) { //without +1
            curr.second = max(curr.second, p.second);
        } else {
            superTranscriptIntervals.push_back(curr);
            curr = p;
        };
    };
    superTranscriptIntervals.push_back(curr);//record last element

    //superTranscript sequences
    uint64 maxSTlen=0;
    for(const auto& p: superTranscriptIntervals) {
        sequenceOfSuperTranscripts.push_back(vector<uint8>(sequenceOfCondensedGenome.begin()+p.first, sequenceOfCondensedGenome.begin()+p.second));
        maxSTlen=max(maxSTlen,p.second-p.first);
    };
    P.inOut->logMain << "Number of superTranscripts = " << superTranscriptIntervals.size() <<";   max length = " <<maxSTlen <<endl;
    
    normalTranscriptSuperTindex.resize(transcriptID.size());
    //for each normal transcript, find superTranscript it belongs to
    uint64 ist=0;
    for (auto &ee: exonLoci) {
        if (ee[exS]>superTranscriptIntervals[ist].second)
            ++ist;
        normalTranscriptSuperTindex[ee[exT]]=ist;
    };
    
    
    normalTranscriptIntervalsInST.resize(transcriptID.size());
    for(uint i = 0; i < transcriptIntervals.size(); i++) {
        normalTranscriptIntervalsInST[i].first = transcriptIntervals[i].first - superTranscriptIntervals[normalTranscriptSuperTindex[i]].first;
        normalTranscriptIntervalsInST[i].second = transcriptIntervals[i].second - superTranscriptIntervals[normalTranscriptSuperTindex[i]].first;
    }
    //////////////////Normal transcripts
    sort(exonLoci.begin(), exonLoci.end(),
         [](const array<uint64, exL>& e1, const array<uint64, exL>& e2) {
             return (e1[exT] < e2[exT]) || ((e1[exT] == e2[exT]) && (e1[exS] < e2[exS]));
         });
    
    sequenceOfNormalTranscripts.resize(transcriptID.size());
    for(uint64 i = 0; i < exonLoci.size(); i++) {
        sequenceOfNormalTranscripts[exonLoci[i][exT]].insert(sequenceOfNormalTranscripts[exonLoci[i][exT]].end(),
                                                             sequenceOfCondensedGenome.begin()+exonLoci[i][exS], sequenceOfCondensedGenome.begin()+exonLoci[i][exE] + 1);
    };

    //output normal transcript sequences
    vector<char> numToCharConverter{'A','C','G','T','N','M'};
    ofstream & tranIdSequences = ofstrOpen(P.pGe.gDir+"/transcriptSequences.fasta",ERROR_OUT, P);
    for(uint64 i = 0; i < sequenceOfNormalTranscripts.size(); i++) {
        tranIdSequences << ">" << transcriptID[i] << "\n";
        for(uint64 j = 0; j < sequenceOfNormalTranscripts[i].size(); j++) {
            tranIdSequences << numToCharConverter[sequenceOfNormalTranscripts[i][j]];
        }
        tranIdSequences << "\n";
    }
    tranIdSequences.close();
    
    //splice junctions, in superTranscript coordinates. Here, sjS is the last base of donor exon, and sjE is the last base of acceptor exon
    for(uint64 i = 1; i < exonLoci.size(); i++) {//exonLoci are still sorted by transcript index and start coordinate
        if (exonLoci[i][exT]==exonLoci[i-1][exT]) {//same transcript
            if (exonLoci[i][exS] > exonLoci[i-1][exE]+1) {//gap >0, otherwise we do not need to record the junction (may need it in the future for other purposes though)
                uint64 sti=normalTranscriptSuperTindex[exonLoci[i][exT]];//ST index
                uint64 sts=superTranscriptIntervals[sti].first; //superTranscript start
                spliceJunctions.emplace_back();//add one element
                spliceJunctions.back()[sjS] = (uint32)(exonLoci[i-1][exE]-sts);
                spliceJunctions.back()[sjE] = (uint32)(exonLoci[i]  [exS]-sts);
                spliceJunctions.back()[sjT] = exonLoci[i][exT];
                spliceJunctions.back()[sjSu] = sti;
            };
        };
    };
    //sort junctions by coordinate 
    sort(spliceJunctions.begin(), spliceJunctions.end(),
     [](const array<uint32, sjL>& e1, const array<uint32, sjL>& e2) {
         return (e1[sjSu] < e2[sjSu]) || ((e1[sjSu] == e2[sjSu]) && (e1[sjE] < e2[sjE]));
     });
    P.inOut->logMain << "Number of splice junctions in superTranscripts = " << spliceJunctions.size() <<endl;

    
    //replace the Genome sequence and all necessary parameters for downstream genome generation
    if (P.pGe.gType=="Transcriptome") {
        genome.concatenateChromosomes(sequenceOfNormalTranscripts, transcriptID, genome.genomeChrBinNbases);
        gtfYes=false; //annotations were converted into transcript sequences and are no longer used
        P.sjdbInsert.yes=false; //actually, it might pe possible to add splice junctions on top of the transcriptome
    } else if (P.pGe.gType=="SuperTranscriptome") {
        vector<string> superTranscriptID(superTranscriptIntervals.size());
        for (uint64 ii=0; ii<superTranscriptID.size(); ii++)
            superTranscriptID[ii]="st" + to_string(ii);
        
        genome.concatenateChromosomes(sequenceOfSuperTranscripts, superTranscriptID, genome.genomeChrBinNbases);
        
        transcriptStrand.resize(superTranscriptIntervals.size(),1);
        
        //sort but exon start position
        sort(exonLoci.begin(), exonLoci.end(),
            [](const array<uint64, exL>& e1, const array<uint64, exL>& e2)
            {
                return e1[exS] < e2[exS];
            }
        );
        //shift exonLoci according to chrStart
        uint64 ist=0;
        for (auto &ee: exonLoci) {
            if (ee[exS]>superTranscriptIntervals[ist].second)
                ++ist;
            ee[exS]+=genome.chrStart[ist]-superTranscriptIntervals[ist].first;
            ee[exE]+=genome.chrStart[ist]-superTranscriptIntervals[ist].first;
        };
    };
}

void Genome::concatenateChromosomes(const vector<vector<uint8>> &vecSeq, const vector<string> &vecName, const uint64 padBin)
{//concatenate chromosome sequence into linear genome, with padding
 //changes genome: chrStart, chrLen, chrName, G  
    nChrReal=vecSeq.size();
    chrLength.resize(nChrReal,0);
    chrStart.resize(nChrReal+1,0);
    chrName=vecName;
    chrNameIndex.clear();
    for (uint32 ii=0; ii<nChrReal; ii++) {
        chrLength[ii]=vecSeq[ii].size();
        chrStart[ii+1]=chrStart[ii]+((chrLength[ii]+1)/padBin+1)*padBin;//+1 makes sure that there is at least one spacer base between chromosomes
        chrNameIndex[chrName[ii]]=ii;
    };

    nGenome=chrStart.back();
    //re-allocate G and copy the sequences
    delete[] G1;
    genomeSequenceAllocate();
    for (uint32 ii=0; ii<nChrReal; ii++) {
        memcpy(G+chrStart[ii],vecSeq[ii].data(),vecSeq[ii].size());
    };
    
    for (uint ii=0;ii<nGenome;ii++) {//- strand
        G[2*nGenome-1-ii]=G[ii]<4 ? 3-G[ii] : G[ii];
    };    
};
