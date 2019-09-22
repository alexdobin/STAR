#include "SuperTranscript.h"
#include "streamFuns.h"

void SuperTranscript::sjCollapse()
{//sort junctions by superTranscript, and then by the end coordinate, and then by the start coordinate
    sort(sj.begin(), sj.end(),
     [](const sjInfo &sj1, const sjInfo &sj2) {
         return ( sj1.super  < sj2.super) || 
                ( sj1.super == sj2.super && sj1.end  < sj2.end ) ||
                ( sj1.super == sj2.super && sj1.end == sj2.end && sj1.start < sj2.start) ;
     });

    //collapse junctions in each superTr, those junctions that belong to different transcript - remove the transcript info
    sjC.resize(seq.size());
    for(uint64 i = 0; i < sj.size(); i++) {
        if (i==0 || sj[i].start!=sj[i-1].start || sj[i].end!=sj[i-1].end || sj[i].super!=sj[i-1].super) //new junction
            sjC[ sj[i].super ].push_back({sj[i].start, sj[i].end});
    };
    
    ofstream & superTrSJ = ofstrOpen(P.pGe.gDir+"/superTranscriptSJcollapsed.tsv", ERROR_OUT, P);
    for(uint64 i = 0; i < sjC.size(); i++) {
        for(auto &sj1 : sjC[i])
            superTrSJ << i <<"\t"<< sj1[0] <<"\t"<< sj1[1] << "\n";
    };
    superTrSJ.close(); 

    P.inOut->logMain << "Number of splice junctions in superTranscripts = " << sj.size() <<endl;
    P.inOut->logMain << "Number of collapsed splice junctions in superTranscripts = " << sjC.size() <<endl;
};

void SuperTranscript::load(char *G, vector<uint64> &chrStart)
{
    seqp.resize(length.size());
    for (uint64 ii=0; ii<length.size(); ii++)
        seqp[ii]=(uint8*)G+chrStart[ii];
    
    ifstream & superTrSJ = ifstrOpen(P.pGe.gDir+"/superTranscriptSJcollapsed.tsv", ERROR_OUT, "SOLUTION: re-generate the genome.", P);
    sjC.resize(length.size());
    uint32 sutr=0,sutr1=0;
    vector<array<uint32,2>> sj1;
    while(superTrSJ >> sutr) {
        if (sutr!=sutr1) {
            sjC[sutr1]=sj1;
            sj1.clear();
            sutr1=sutr;
        };
        sj1.emplace_back();
        superTrSJ >> sj1.back()[0] >> sj1.back()[1];
    };
    superTrSJ.close();
    
};
