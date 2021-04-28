#include "SoloReadFeature.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ReadAnnotations.h"
#include "SoloReadBarcode.h"

class ReadSoloFeatures {
public:
    uint32 gene;
    vector<uint32> geneMult;
    vector<array<uint64,2>> sj;
    bool sjAnnot;
    uint32 indAnnotTr; //index of the annotated transcript
    Transcript **alignOut;
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot);

void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (pSolo.type==0 || soloBar.cbMatch<0)
        return;
       
    ReadSoloFeatures reFe;
    reFe.alignOut=alignOut;
    
    uint32 nFeat=0; //number of features in this read (could be >1 for SJs)
    if (nTr==0) {//unmapped
        stats.V[stats.nUnmapped]++;
        
    } else {
        switch (featureType) {
            case SoloFeatureTypes::Gene :
            case SoloFeatureTypes::GeneFull : 
                {
                    set<uint32> *readGe=( featureType==SoloFeatureTypes::Gene 
                                                      ? &readAnnot.geneConcordant : &readAnnot.geneFull ); //for Gene or GeneFull
                    if (readGe->size()==0) {//check genes
                        stats.V[stats.nNoFeature]++;//no gene
                    } else if (readGe->size()>1) {
                        stats.V[stats.nAmbigFeature]++;//multigene
                        if (nTr>1)
                            stats.V[stats.nAmbigFeatureMultimap]++;//multigene caused by multimapper
                            
                        if (pSolo.multiMap.yes.multi) {//output multimappers
                            reFe.geneMult.resize(readGe->size());
                            uint32 ii=0;
                            for (auto &g: *readGe) {//set high bit to mark multimappers
                                reFe.geneMult[ii] = g | geneMultMark;
                                ++ii;
                            };
                            
                            reFe.indAnnotTr = ( featureType==SoloFeatureTypes::Gene ? readAnnot.geneConcordantTr : readAnnot.geneFullTr );//TODO deal carefully with start/end collapsing for SmartSeq                  
                            nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot);
                        };
                    } else {//good gene
                        reFe.gene = *readGe->begin();
                        reFe.indAnnotTr = ( featureType==SoloFeatureTypes::Gene ? readAnnot.geneConcordantTr : readAnnot.geneFullTr );                        
                        nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (uint64)-1), featureType, soloBar, reFe, readAnnot);
                    };
                };
                break;
        
            case SoloFeatureTypes::SJ : 
                if (nTr>1) {//reject all multimapping junctions
                    stats.V[stats.nAmbigFeatureMultimap]++;
                } else if (readAnnot.geneConcordant.size()>1){//for SJs, still check genes, no feature if multi-gene
                    stats.V[stats.nAmbigFeature]++;
                } else {//one gene or no gene
                    alignOut[0]->extractSpliceJunctions(reFe.sj, reFe.sjAnnot);
                    if ( reFe.sj.empty() || (reFe.sjAnnot && readAnnot.geneConcordant.size()==0) ) {//no junctions, or annotated junction but no gene (i.e. read does not fully match transcript)
                        stats.V[stats.nNoFeature]++;
                    } else {//good junction
                        reFe.indAnnotTr=0;                        
                        nFeat = outputReadCB(streamReads, (readIndexYes ? iRead : (uint64)-1), featureType, soloBar, reFe, readAnnot);
                    };
                };                  
                break;
        
            case SoloFeatureTypes::Transcript3p : 
                if (readAnnot.transcriptConcordant.size()==0 || soloBar.cbMatch>1) {//do not record ambiguous CB  
                    stats.V[stats.nNoFeature]++;
                } else {
                    nFeat = outputReadCB(streamReads, iRead, featureType, soloBar, reFe, readAnnot);
                };                
                if (readAnnot.transcriptConcordant.size()==1 && readAnnot.transcriptConcordant[0][1] < transcriptDistCount.size()) {
                    //read maps to one transcript - use for distTTS distribution function
                    transcriptDistCount[readAnnot.transcriptConcordant[0][1]]++;
                };
                
                break;

            case SoloFeatureTypes::Velocyto :
                //different record: iRead, nTr, tr1, type1, tr2, type2 ...
                if (readAnnot.trVelocytoType.size()>0) {//otherwise, no gene
                    
                    sort(readAnnot.trVelocytoType.begin(), readAnnot.trVelocytoType.end(),
                         [](const trTypeStruct &t1, const trTypeStruct &t2) {return t1.tr < t2.tr;});

                    *streamReads << iRead <<' '<< readAnnot.trVelocytoType.size();
                    for (auto &tt: readAnnot.trVelocytoType)
                         *streamReads <<' '<< tt.tr <<' '<< (uint32) tt.type;
                    *streamReads <<'\n';
                    nFeat=1;
                } else {
                    stats.V[stats.nNoFeature]++;
                };
                break; //no need to go with downstream processing                
                
        };//switch (featureType)
    };//if (nTr==0)
    
    if (nFeat==0 && readInfoYes) {//no feature, but readInfo requested
        outputReadCB(streamReads, iRead, (uint32)-1, soloBar, reFe, readAnnot);
    };
    
    if (nFeat==0)
        return; //no need to record the number of reads per CB
    
    if (pSolo.cbWLyes) {//WL
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nFeat;
    } else {//no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nFeat;
    };
    
    return;
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, SoloReadBarcode &soloBar, const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot)
{   
    /*format of the temp output file
     * UMI [iRead] type feature* cbMatchString
     *             0=exact match, 1=one non-exact match, 2=multipe non-exact matches
     *                   gene or sj[0] sj[1]
     *                         CB or nCB {CB Qual, ...}
     */
    
    if (soloBar.pSolo.type==soloBar.pSolo.SoloTypes::SmartSeq && featureType!=-1) {//need to calculate "UMI" from align start/end
        soloBar.umiB=reFe.alignOut[reFe.indAnnotTr]->chrStartLengthExtended();
    };
    
    uint64 nout=1;
    
    switch (featureType) {
        case -1 :
            //no feature, output for readInfo
            *streamOut << soloBar.umiB <<' '<< iRead <<' '<< -1 <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
            break;
            
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
            
            if (reFe.geneMult.size()==0) {
                //just gene id
                *streamOut << soloBar.umiB <<' ';//UMI
                if ( iRead != (uint64)-1 )
                    *streamOut << iRead <<' ';//iRead
                *streamOut << reFe.gene <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';;
            } else {
                for (auto &g : reFe.geneMult) {
                    *streamOut << soloBar.umiB <<' ';//UMI
                    *streamOut << iRead <<' ';//iRead is always output for multiGene
                    *streamOut << g <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';;
                };
                nout = reFe.geneMult.size();
            };
            break;
            
        case SoloFeatureTypes::SJ :
            //sj - two numbers, multiple sjs per read
            for (auto &sj : reFe.sj) {
                *streamOut << soloBar.umiB <<' ';//UMI
                if ( iRead != (uint64)-1 )
                    *streamOut << iRead <<' ';//iRead            
                *streamOut << sj[0] <<' '<< sj[1] <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
            };
            nout=reFe.sj.size();
            break;
            
        case SoloFeatureTypes::Transcript3p :
            //transcript,distToTTS structure    
            *streamOut << soloBar.cbMatchString <<' ';            
            *streamOut << soloBar.umiB <<' ';
            *streamOut << readAnnot.transcriptConcordant.size();
            for (auto &tt: readAnnot.transcriptConcordant) {
                *streamOut <<' '<< tt[0] <<' '<< tt[1];
            };
            if ( iRead != (uint64)-1 )
                *streamOut  <<' '<< iRead;//iRead
            *streamOut  <<'\n';
            nout=1;

            break;
    }; //switch (featureType)
    
    return nout;
};
