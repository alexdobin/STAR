#include "samHeaders.h"
#include <iostream>
#include "BAMfunctions.h"

void samHeaders(Parameters &P, Genome &genomeOut, Transcriptome &transcriptomeMain) 
{
    /////////////////////////////////////////////////////////////////////////////////// transcriptome BAM header
    if ( P.quant.trSAM.bamYes ) {//header for transcriptome BAM
        ostringstream samHeaderStream;
        vector <uint> trlength;
        for (uint32 ii=0;ii<transcriptomeMain.trID.size();ii++) {
            uint32 iex1=transcriptomeMain.trExI[ii]+transcriptomeMain.trExN[ii]-1; //last exon of the transcript
            trlength.push_back(transcriptomeMain.exLenCum[iex1]+transcriptomeMain.exSE[2*iex1+1]-transcriptomeMain.exSE[2*iex1]+1);
            samHeaderStream << "@SQ\tSN:"<< transcriptomeMain.trID.at(ii) <<"\tLN:"<<trlength.back()<<"\n";
        };
        for (uint32 ii=0;ii<P.outSAMattrRGlineSplit.size();ii++) {//@RG lines
            samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) <<"\n";
        };
        outBAMwriteHeader(P.inOut->outQuantBAMfile,samHeaderStream.str(),transcriptomeMain.trID,trlength);
    };
    
    //////////////////////////////////////////////////////////////////////////////// main headers
    if (P.outSAMmode == "None" || P.outSAMtype[0] == "None") {//no SAM output
        return;
    };
    
    ostringstream samHeaderStream;
    
    for (uint ii=0;ii<genomeOut.nChrReal;ii++) {
        samHeaderStream << "@SQ\tSN:"<< genomeOut.chrName.at(ii) <<"\tLN:"<<genomeOut.chrLength[ii]<<"\n";
    };

    genomeOut.chrNameAll=genomeOut.chrName;
    genomeOut.chrLengthAll=genomeOut.chrLength;
    {//add exra references
        ifstream extrastream (P.pGe.gDir + "/extraReferences.txt");
        while (extrastream.good()) {
            string line1;
            getline(extrastream,line1);
            istringstream stream1 (line1);
            string field1;
            stream1 >> field1;//should check for @SQ

            if (field1!="") {//skip blank lines
                samHeaderStream << line1 <<"\n";

                stream1 >> field1;
                genomeOut.chrNameAll.push_back(field1.substr(3));
                stream1 >> field1;
                genomeOut.chrLengthAll.push_back((uint) stoll(field1.substr(3)));
            };
        };
        extrastream.close();
    };

    if (P.outSAMheaderPG.at(0)!="-") {
        samHeaderStream << P.outSAMheaderPG.at(0);
        for (uint ii=1;ii<P.outSAMheaderPG.size(); ii++) {
            samHeaderStream << "\t" << P.outSAMheaderPG.at(ii);
        };
        samHeaderStream << "\n";
    };

    samHeaderStream << "@PG\tID:STAR\tPN:STAR\tVN:" << STAR_VERSION <<"\tCL:" << P.commandLineFull <<"\n";

    if (P.outSAMheaderCommentFile!="-") {
        ifstream comstream (P.outSAMheaderCommentFile);
        while (comstream.good()) {
            string line1;
            getline(comstream,line1);
            if (line1.find_first_not_of(" \t\n\v\f\r")!=std::string::npos) {//skip blank lines
                samHeaderStream << line1 <<"\n";
            };
        };
        comstream.close();
    };


    for (uint32 ii=0;ii<P.outSAMattrRGlineSplit.size();ii++) {//@RG lines
        samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) <<"\n";
    };


    samHeaderStream <<  "@CO\t" <<"user command line: " << P.commandLine <<"\n";

    samHeaderStream << P.samHeaderExtra;

    if (P.outSAMheaderHD.at(0)!="-") {
        P.samHeaderHD = P.outSAMheaderHD.at(0);
        for (uint ii=1;ii<P.outSAMheaderHD.size(); ii++) {
            P.samHeaderHD +="\t" + P.outSAMheaderHD.at(ii);
        };
    } else {
        P.samHeaderHD = "@HD\tVN:1.4";
    };


    P.samHeader=P.samHeaderHD+"\n"+samHeaderStream.str();
    //for the sorted BAM, need to add SO:cooridnate to the header line
    P.samHeaderSortedCoord=P.samHeaderHD + (P.outSAMheaderHD.size()==0 ? "" : "\tSO:coordinate") + "\n" + samHeaderStream.str();

    if (P.outSAMbool) {//
        *P.inOut->outSAM << P.samHeader;
    };
    if (P.outBAMunsorted){
        outBAMwriteHeader(P.inOut->outBAMfileUnsorted,P.samHeader,genomeOut.chrNameAll,genomeOut.chrLengthAll);
    };
};