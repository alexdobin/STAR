#include "SuperTranscriptome.h"
#include "streamFuns.h"

void SuperTranscriptome::sjCollapse()
{//sort junctions by superTranscript, and then by the end coordinate, and then by the start coordinate
    sort(sj.begin(), sj.end(),
     [](const sjInfo &sj1, const sjInfo &sj2) {
         return ( sj1.super  < sj2.super) || 
                ( sj1.super == sj2.super && sj1.start < sj2.start ) ||
                ( sj1.super == sj2.super && sj1.start == sj2.start && sj1.end < sj2.end ) ;
     });

    //collapse junctions in each superTr, those junctions that belong to different transcript - remove the transcript info
    vector<vector<array<uint32,2>>> sjCollapsed;
    sjCollapsed.resize(seq.size());
    for(uint32 i = 0; i < sj.size(); i++) {
        if (i==0 || sj[i].start!=sj[i-1].start || sj[i].end!=sj[i-1].end || sj[i].super!=sj[i-1].super) //new junction
            sjCollapsed[ sj[i].super ].push_back({sj[i].start, sj[i].end});
    };
    
    ofstream & superTrSJstream = ofstrOpen(P.pGe.gDir+"/superTranscriptSJcollapsed.tsv", ERROR_OUT, P);
    for(uint64 i = 0; i < sjCollapsed.size(); i++) {
        for(auto &sj1 : sjCollapsed[i])
            superTrSJstream << i <<"\t"<< sj1[0] <<"\t"<< sj1[1] << "\n";
    };
    superTrSJstream.close();

    P.inOut->logMain << "Number of splice junctions in superTranscripts = " << sj.size() <<endl;
    P.inOut->logMain << "Number of collapsed splice junctions in superTranscripts = " << sjCollapsed.size() <<endl;
};

void SuperTranscriptome::load(char *G, vector<uint64> &chrStart, vector<uint64> &chrLength)
{//load superTranscript seqs and
    N=chrLength.size();
    superTrs.resize(N);
    for (uint32 ii=0; ii<N; ii++) {
        superTrs[ii].length=chrLength[ii];
        superTrs[ii].seqP=(uint8*)G+chrStart[ii];
    };
    
    ifstream & superTrSJstream = ifstrOpen(P.pGe.gDir+"/superTranscriptSJcollapsed.tsv", ERROR_OUT, "SOLUTION: re-generate the genome.", P);
    
    uint32 sutr=0,sutr1=0;
    vector<array<uint32,3>> sjC1;
    sjNmax=0;
    sjDonorNmax=0;
    vector<uint32> sjDonor1;
    
    while(true) {//load junctions, assume they are sorted by donor coordinate
        
        superTrSJstream >> sutr;
        bool inGood = superTrSJstream.good();
        
        if (sutr!=sutr1 || !inGood) {//new suTr
            //sort sj1 by acceptor position
            sort(sjC1.begin(), sjC1.end(),
                [](const array<uint32,3> &sj1, const array<uint32,3> &sj2) {
                 return ( sj1[1] < sj2[1] ) ||
                        ( sj1[1] == sj2[1]  && sj1[0] < sj2[0] );
                });

            superTrs[sutr1].sjC=sjC1;
            superTrs[sutr1].sjDonor=sjDonor1;
            
            if (sjNmax < sjC1.size())
                sjNmax=sjC1.size();
            if (sjDonorNmax < sjDonor1.size()) {
                sjDonorNmax=sjDonor1.size();
            };
            
            sjC1.clear();
            sjDonor1.clear();
            sutr1=sutr;
            
            if (!inGood) {
                break; //end of file
            };
        };
        
        uint32 sjd, sja;
        superTrSJstream >> sjd >> sja;
        
        if (sjDonor1.empty() || sjDonor1.back() < sjd) 
            sjDonor1.push_back(sjd);
        
        sjC1.push_back({sjd,sja,(uint32)sjDonor1.size()-1});//record donor, acceptor, and position of the donor in sjDonor
    };
    superTrSJstream.close();
    
    P.inOut->logMain << "Max number of splice junctions in a superTranscript = " << sjNmax <<endl;
    P.inOut->logMain << "Max number of donor sites in a superTranscript = " << sjDonorNmax <<endl;
};
