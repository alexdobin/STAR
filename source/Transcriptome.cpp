#include "Transcriptome.h"
#include "streamFuns.h"
#include "GlobalVariables.h"

Transcriptome::Transcriptome (Parameters &Pin) : P(Pin){

    trInfoDir = P.pGe.sjdbGTFfile=="-" ? P.pGe.gDir : P.sjdbInsert.outDir; //if GTF file is given at the mapping stage, it's always used for transcript info

    if ( P.quant.trSAM.yes ) {//load exon-transcript structures
        //load tr and ex info

        ifstream & trinfo = ifstrOpen(trInfoDir+"/transcriptInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step",P);

        trinfo >> nTr;
        trS=new uint [nTr];
        trE=new uint [nTr];
        trEmax=new uint [nTr];
        trExI=new uint32 [nTr];
        trExN=new uint16 [nTr];
        trStr=new uint8 [nTr];
        trID.resize(nTr);

        for (uint32 itr=0; itr<nTr; itr++) {
            uint16 str1;
            trinfo >> trID[itr] >> trS[itr] >> trE[itr] >> trEmax[itr] >> str1 >> trExN[itr] >> trExI[itr];
            trStr[itr]=str1;
        };
        P.inOut->logMain << "Loaded transcript database, nTr="<<nTr<<endl;

        trinfo.close();

        ifstream & exinfo = ifstrOpen(trInfoDir+"/exonInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step", P);

        exinfo >> nEx;
        exSE = new uint32 [2*nEx];
        exLenCum = new uint32 [nEx];
        for (uint32 iex=0; iex<nEx; iex++) {
            exinfo >> exSE[2*iex] >> exSE[2*iex+1] >> exLenCum[iex]; //reading all elements one after another
        };
        P.inOut->logMain << "Loaded exon database, nEx="<<nEx<<endl;

        exinfo.close();

        //
    };

    if ( P.quant.geCount.yes ) {//load exon-gene structures
        ifstream & exinfo = ifstrOpen(trInfoDir+"/exonGeTrInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step", P);
        exinfo >> exG.nEx;

        exG.s=new uint64[exG.nEx];
        exG.e=new uint64[exG.nEx];
        exG.eMax=new uint64[exG.nEx];
        exG.str=new uint8[exG.nEx];
        exG.g=new uint32[exG.nEx];
        exG.t=new uint32[exG.nEx];


        for (uint ii=0;ii<exG.nEx;ii++) {
            int str1;
            exinfo >> exG.s[ii] >> exG.e[ii] >> str1 >> exG.g[ii] >> exG.t[ii];
            exG.str[ii] = (uint8) str1;
        };
        exinfo.close();

        //calculate eMax
        exG.eMax[0]=exG.e[0];
        for (uint iex=1;iex<exG.nEx;iex++) {
            exG.eMax[iex]=max(exG.eMax[iex-1],exG.e[iex]);
        };

        ifstream & geStream = ifstrOpen(trInfoDir+"/geneInfo.tab", ERROR_OUT, "SOLUTION: utilize --sjdbGTFfile /path/to/annotantions.gtf option at the genome generation step or mapping step", P);
        geStream >> nGe;
        geID.resize(nGe);
        for (uint ii=0;ii<nGe;ii++) {
            geStream >> geID[ii];
        };
        geStream.close();


    };
};

void Transcriptome::quantsAllocate() {
    if ( P.quant.geCount.yes ) {
        quants = new Quantifications (nGe);
    };
};

void Transcriptome::quantsOutput() {
    ofstream qOut(P.quant.geCount.outFile);

    qOut << "N_unmapped";
    for (int itype=0; itype<quants->geneCounts.nType; itype++)
    {
        qOut << "\t" <<g_statsAll.unmappedAll;
    };
    qOut << "\n";

    qOut << "N_multimapping";
    for (int itype=0; itype<quants->geneCounts.nType; itype++)
    {
        qOut << "\t" <<quants->geneCounts.cMulti ;
    };
    qOut << "\n";

    qOut << "N_noFeature";
    for (int itype=0; itype<quants->geneCounts.nType; itype++)
    {
        qOut << "\t" <<quants->geneCounts.cNone[itype];
    };
    qOut << "\n";

    qOut << "N_ambiguous";
    for (int itype=0; itype<quants->geneCounts.nType; itype++)
    {
        qOut << "\t" <<quants->geneCounts.cAmbig[itype];
    };
    qOut << "\n";

    for (uint32 ig=0; ig<nGe; ig++)
    {
        qOut << geID[ig];
        for (int itype=0; itype<quants->geneCounts.nType; itype++)
        {
            qOut << "\t" <<quants->geneCounts.gCount[itype][ig];
        };
        qOut << "\n";
    };
    qOut.close();
};
