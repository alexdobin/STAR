#include "SuperTranscript.h"
#include "streamFuns.h"

void SuperTranscript::sjCollapse()
{//sort junctions by superTranscript, and then by the end coordinate, and then by the start coordinate
    sort(sj.begin(), sj.end(),
     [](const sjInfo &sj1, const sjInfo &sj2) {
         return ( sj1.super  < sj2.super) || 
                ( sj1.super == sj2.super && sj1.start < sj2.start ) ||
                ( sj1.super == sj2.super && sj1.start == sj2.start && sj1.end < sj2.end ) ;
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
{//load superTranscript seqs and
    seqp.resize(N);
    for (uint64 ii=0; ii<N; ii++)
        seqp[ii]=(uint8*)G+chrStart[ii];
    
    ifstream & superTrSJ = ifstrOpen(P.pGe.gDir+"/superTranscriptSJcollapsed.tsv", ERROR_OUT, "SOLUTION: re-generate the genome.", P);
    
    sjC.resize(N);
    sjDonor.resize(N);
    
    uint32 sutr=0,sutr1=0;
    vector<array<uint32,3>> sjC1;
    sjNmax=0;
    vector<uint32> sjDonor1;
    
    while(superTrSJ >> sutr) {//load junctions, assume they are sorted by donor coordinate
        if (sutr!=sutr1) {//new suTr
            //sort sj1 by acceptor position
            sort(sjC1.begin(), sjC1.end(),
                [](const array<uint32,3> &sj1, const array<uint32,3> &sj2) {
                 return ( sj1[1] < sj2[1] ) ||
                        ( sj1[1] == sj2[1]  && sj1[0] < sj2[0] );
                });

            sjC[sutr1]=sjC1;
            
            if (sjNmax < sjC1.size())
                sjNmax=sjC1.size();
            
            sjC1.clear();
            sutr1=sutr;
        };
        
        uint32 sjd, sja;
        superTrSJ >> sjd >> sja;
        
        if (sjDonor1.back() < sjd) 
            sjDonor1.push_back(sjd);
        
        sjC1.push_back({sjd,sja,(uint32)sjDonor1.size()-1});//record donor, acceptor, and position of the donor in sjDonor
    };
    superTrSJ.close();
};
