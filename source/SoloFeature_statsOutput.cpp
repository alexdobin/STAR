#include "SoloFeature.h"
#include "streamFuns.h"
#include "Stats.h"
#include "GlobalVariables.h"

void SoloFeature::statsOutput()
{
    ofstream &strOut=ofstrOpen(outputPrefix+"Summary.csv", ERROR_OUT, P);
    //Sequencing
    strOut << "Number of Reads," << g_statsAll.readN <<'\n';
    strOut << "Reads With Valid Barcodes," << 1.0 - double( readBarSum->stats.numInvalidBarcodes() + readFeatSum->stats.numInvalidBarcodes() )/g_statsAll.readN <<'\n';
    strOut << "Sequencing Saturation," << readFeatSum->stats.numSequencingSaturation() <<'\n';
    
    if (pSolo.type != pSolo.SoloTypes::SmartSeq) {//qualHist for CB+UMI
        uint64 q30=0, ntot=0;
        for (uint32 ix=0; ix<256; ix++) {
            ntot += readBarSum->qualHist[ix];
            if (ix >= (P.readQualityScoreBase + 30))
                q30 += readBarSum->qualHist[ix];
        };
           
        strOut << "Q30 Bases in CB+UMI,"   << double(q30)/ntot <<'\n';
    };
    
    {//qualHist for RNA
        uint64 q30=0, ntot=0;
        for (int ichunk=0; ichunk<P.runThreadN; ichunk++) {
            for (uint32 imate=0; imate<P.readNmates; imate++) {
                for (uint32 ix=0; ix<256; ix++) {
                    ntot += RAchunk[ichunk]->RA->qualHist[imate][ix];
                    if (ix >= (P.readQualityScoreBase + 30))
                        q30 += RAchunk[ichunk]->RA->qualHist[imate][ix];
                };
            };
        };
        strOut << "Q30 Bases in RNA read," << double(q30)/ntot <<'\n';
    };    
    
    strOut << "Reads Mapped to Genome: Unique+Multiple," << double(g_statsAll.mappedReadsU+g_statsAll.mappedReadsM)/g_statsAll.readN <<'\n';
    strOut << "Reads Mapped to Genome: Unique," << double(g_statsAll.mappedReadsU)/g_statsAll.readN <<'\n';
    
    string mapfeat=SoloFeatureTypes::Names[featureType];
    
    strOut << "Reads Mapped to "<< mapfeat << ": Unique+Multiple " << SoloFeatureTypes::Names[featureType] <<",";
    if (pSolo.multiMap.yes.multi) {
        strOut << double( readFeatSum->stats.numMappedToTranscriptome() )/g_statsAll.readN <<'\n';
    } else {
        strOut << "NoMulti\n";
    };

    strOut << "Reads Mapped to "<< mapfeat << ": Unique " << SoloFeatureTypes::Names[featureType] <<"," << double( readFeatSum->stats.numMappedToTranscriptomeUnique() )/g_statsAll.readN <<'\n';
    
    if (pSolo.cellFilter.type[0]!="None"
        && (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull
         || featureType==SoloFeatureTypes::GeneFull_ExonOverIntron || featureType==SoloFeatureTypes::GeneFull_Ex50pAS)) {
        //if (pSolo.cellFilter.type[0]=="CellRanger2.2") 
        {
            strOut << "Estimated Number of Cells," << filteredCells.nCells <<'\n';

            strOut << "Unique Reads in Cells Mapped to " << SoloFeatureTypes::Names[featureType] << "," << filteredCells.nReadInCellsUnique <<'\n';
            strOut << "Fraction of Unique Reads in Cells," << double(filteredCells.nReadInCellsUnique) / readFeatSum->stats.numMappedToTranscriptomeUnique() <<'\n';
            strOut << "Mean Reads per Cell," << filteredCells.meanReadPerCellUnique <<'\n';
            strOut << "Median Reads per Cell," << filteredCells.medianReadPerCellUnique <<'\n';

            strOut << "UMIs in Cells," << filteredCells.nUMIinCells <<'\n';
            strOut << "Mean UMI per Cell," << filteredCells.meanUMIperCell <<'\n';
            strOut << "Median UMI per Cell," << filteredCells.medianUMIperCell <<'\n';    

            strOut << "Mean "   << SoloFeatureTypes::Names[featureType] << " per Cell," << filteredCells.meanGenePerCell <<'\n';
            strOut << "Median " << SoloFeatureTypes::Names[featureType] << " per Cell," << filteredCells.medianGenePerCell <<'\n';    
            strOut << "Total "  << SoloFeatureTypes::Names[featureType] << " Detected," << filteredCells.nGeneDetected <<'\n';    
        };

        //output UMI per cell, sorted
        ofstream &strOutUMIperCell = ofstrOpen(outputPrefix+"UMIperCellSorted.txt", ERROR_OUT, P);

        for (auto & n : nUMIperCBsorted) {
            if (n==0)
                break;
            strOutUMIperCell << n <<'\n';
        };
        strOutUMIperCell.close();
    };
    
    strOut.close();

    ///////////////////////////////////////////////// readStatsOutput
    if (pSolo.readStatsYes[featureType]) {
        ofstream &strOut=ofstrOpen(outputPrefix+"CellReads.stats", ERROR_OUT, P);

        strOut << "CB";
        for (auto &sn: readFlagCounts.statNames)
            strOut <<"\t"<< sn;
        strOut <<"\t"<< "nUMIunique" <<"\t"<< "nGenesUnique" <<"\t"<< "nUMImulti" <<"\t"<< "nGenesMulti";
        strOut << '\n';

        strOut << "CBnotInPasslist";
        for (auto &cc: readFlagCounts.flagCountsNoCB)
            strOut <<"\t"<< cc;
        strOut <<"\t0\t0\t0\t0\n";

        for (auto &cbc: readFlagCounts.flagCounts) {
            strOut << pSolo.cbWLstr[cbc.first];

            for (uint32 ib=0; ib<readFlagCounts.nBits; ib++)
                strOut <<'\t'<< cbc.second[ib];

            if (indCBwl[cbc.first]==(uint32)-1) {//this CB did not have any reads mapped to features
                strOut <<'\t'<< 0 <<'\t'<< 0 <<'\t'<< 0 <<'\t'<< 0;
            } else {
                strOut <<'\t'<< nUMIperCB[indCBwl[cbc.first]] <<'\t'<< nGenePerCB[indCBwl[cbc.first]];
                if (nUMIperCBmulti.size()==0) {//no multi
                    strOut <<'\t'<< 0 <<'\t'<< 0;
                } else {
                    strOut <<'\t'<< nUMIperCBmulti[indCBwl[cbc.first]] <<'\t'<< nGenePerCBmulti[indCBwl[cbc.first]];
                };
            };
            strOut <<'\n';
        };
        strOut.close();
    };

};
