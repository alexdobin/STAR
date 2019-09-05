#include "SoloReadFeature.h"
#include "serviceFuns.cpp"
#include "SequenceFuns.h"
#include "ReadAnnotations.h"
#include "SoloReadBarcode.h"

class ReadSoloFeatures {
public:
    uint32 gene;
    vector<array<uint64,2>> sj;
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, const SoloReadBarcode &soloBar, const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot);

void SoloReadFeature::record(SoloReadBarcode &soloBar, uint nTr, Transcript *alignOut, uint64 iRead, ReadAnnotations &readAnnot)
{
    if (pSolo.type==0 || soloBar.cbMatch<0)
        return;

    ReadSoloFeatures reFe;
    
    set<uint32> *readGe=&readAnnot.geneConcordant; //for featureType==0

    bool readFeatYes=true;
    //calculate feature
    if (nTr==0) {//unmapped
        stats.V[stats.nUnmapped]++;
        readFeatYes=false;
        
    } else if (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull) {//genes
        //check genes, return if no gene of multimapping
        if (featureType==SoloFeatureTypes::GeneFull) {
            readGe = &readAnnot.geneFull;
        };
        if (readGe->size()==0) {
            stats.V[stats.nNoFeature]++;
            readFeatYes=false;
        };
        if (readGe->size()>1) {
            stats.V[stats.nAmbigFeature]++;
            if (nTr>1)
                stats.V[stats.nAmbigFeatureMultimap]++;
            readFeatYes=false;
        };
        if (readFeatYes)
            reFe.gene=*readGe->begin();
        
    } else if (featureType==SoloFeatureTypes::SJ) {//SJs
        if (nTr>1) {//reject all multimapping junctions
            stats.V[stats.nAmbigFeatureMultimap]++;
            readFeatYes=false;
        } else {//for SJs, still check genes, return if multi-gene
            if (readAnnot.geneConcordant.size()>1) {
                stats.V[stats.nAmbigFeature]++;
                readFeatYes=false;
            } else {//one gene or no gene
                bool sjAnnot;
                alignOut->extractSpliceJunctions(reFe.sj, sjAnnot);
                if ( reFe.sj.empty() || (sjAnnot && readAnnot.geneConcordant.size()==0) ) {//no junctions, or annotated junction but no gene (i.e. read does not fully match transcript)
                    stats.V[stats.nNoFeature]++;
                    readFeatYes=false;
                };
            };
        };
        
    } else if (featureType==SoloFeatureTypes::Transcript3p) {//transcripts
        if (readAnnot.transcriptConcordant.size()==0) {
            stats.V[stats.nNoFeature]++;
            readFeatYes=false;
        };
    };
    
    if (!readFeatYes && !readInfoYes) //no feature, and no readInfo requested
        return;
    
    if (!readInfoYes)
        iRead=(uint64)-1;

    uint32 nfeat = outputReadCB(streamReads, iRead, (readFeatYes ? featureType : -1), soloBar, reFe, readAnnot);
    if (pSolo.cbWLsize>0) {//WL
        for (auto &cbi : soloBar.cbMatchInd)
            cbReadCount[cbi] += nfeat;
    } else {//no WL
        cbReadCountMap[soloBar.cbMatchInd[0]] += nfeat;
    };
};

uint32 outputReadCB(fstream *streamOut, const uint64 iRead, const int32 featureType, const SoloReadBarcode &soloBar, const ReadSoloFeatures &reFe, const ReadAnnotations &readAnnot)
{   
    //format of the temp output file
    // UMI [iRead] type feature* cbMatchString
    //             0=exact match, 1=one non-exact match, 2=multipe non-exact matches
    //                   gene or sj[0] sj[1]
    //                           CB or nCB {CB Qual, ...}
    
    uint64 nout=1;
    if (featureType==-1) {//no feature, output for readInfo
        *streamOut << soloBar.umiB <<' '<< iRead <<' '<< -1 <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';;
        
    } else if (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull) {//Gene, GeneFull
        *streamOut << soloBar.umiB <<' ';//UMI
        if ( iRead != (uint64)-1 )
            *streamOut << iRead <<' ';//iRead
        *streamOut << reFe.gene <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';;
        
    } else if (featureType==SoloFeatureTypes::SJ) {//sjs
        for (auto &sj : reFe.sj) {
            *streamOut << soloBar.umiB <<' ';//UMI
            if ( iRead != (uint64)-1 )
                *streamOut << iRead <<' ';//iRead            
            *streamOut << sj[0] <<' '<< sj[1] <<' '<< soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';;
        };
        nout=reFe.sj.size();
        
    } else if (featureType==SoloFeatureTypes::Transcript3p) {//transcript,distToTTS structure
        *streamOut << soloBar.umiB <<' ';//UMI
        if ( iRead != (uint64)-1 )
            *streamOut << iRead <<' ';//iRead
        *streamOut << readAnnot.transcriptConcordant.size()<<' ';
        for (auto &tt: readAnnot.transcriptConcordant) {
            *streamOut << tt[0] <<' '<< tt[1] <<' ';
        };
        *streamOut << soloBar.cbMatch <<' '<< soloBar.cbMatchString <<'\n';
    };
    
    return nout;
};
