#include "Genome.h"
#include "SequenceFuns.h"
#include "streamFuns.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"

void Genome::transformGenome(vector<array<uint64,exL>> &exonLoci) 
{
    if (pGe.transform.typeString=="None") {
        pGe.transform.type=0;
    } else if (pGe.transform.typeString=="Haploid") {
        pGe.transform.type=1;
    } else if (pGe.transform.typeString=="Diploid") {
        pGe.transform.type=2;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: unrecognized option in --outTransformType=" << pGe.transform.typeString << "\n";
        errOut << "SOLUTION: use one of the allowed values of --outWigType : 'None' or 'Haploid' or 'Diploid' \n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    
    if (pGe.transform.type==0)
        return;
   
//     typedef struct {
//         int32 len;//0: SNV, <0: deletion; >0: insertion
//         bool ref; //true if matches reference
//         string seq;//sequence
//     } AlleleInfo;
    
    typedef struct {
        uint64 pos;
        int32 len;//0: SNV, <0: deletion; >0: insertion
        array<string,2> seq;//sequence for SNV and insertions, empty for deletions
        //array<AlleleInfo,2> info;
    } VariantInfo;
        
    map<string,vector<VariantInfo>> vcfVariants;
    {//load VCF file: per chromosome, 1 or 2 haplotypes
        ifstream &vcfStream = ifstrOpen(pGe.transform.vcfFile, ERROR_OUT, "SOLUTION: check the path and permissions of the VCF file: "+pGe.transform.vcfFile, P);
        string vcfLine;
        while (std::getline (vcfStream, vcfLine)) {
            string chr,id, ref, alt, dummy, sample;
            uint64 pos;
            
            istringstream vcfLineStream(vcfLine);
            vcfLineStream >> chr;
            
            if (chr.at(0)=='#')
                continue;
            if (chrNameIndex.count(chr)==0) {//chr not in Genome
                P.inOut->logMain << "WARNING: while processing varVCFfile file=" << P.var.vcfFile <<": chromosome '"<<chr<<"' not found in Genome fasta file\n";
                continue;
            };
            if (ref=="*" || alt=="*") {
                P.inOut->logMain << "WARNING: VCF: * allele"<<vcfLine<<endl;
                continue;
            };
            
            vcfLineStream >> pos >> id >> ref >> alt >> dummy >> dummy >> dummy >> dummy >> sample;
            //we assume ref can only be one sequence

            vector <string> altV;
            splitString(alt,',',altV);
            
            if (pGe.transform.type==1) {//use first alt allele TODO make warning if there are 2 alts?
                alt=altV[0];
                
                //AlleleInfo alt1={(int32)alt.size()-(int32)ref.size(), false, alt};
                VariantInfo var1={pos, (int32)alt.size()-(int32)ref.size(), {ref,alt}};
                vcfVariants[chr].push_back(var1);
            };
        };
        vcfStream.close();
    };
    //haploid: insert alternative alleles into genome sequence, create conversion-block file
    {   //recalculate chrLength
        vector<uint64> chrStart1(chrStart), chrLength1(chrLength);
        vector<bool> chtTransformYes(chrStart.size(), false);
        for (uint32 ichr=0; ichr<chrLength.size(); ichr++) {
            
            if (vcfVariants.count(chrName[ichr])==0)
                continue;
            
            chtTransformYes[ichr]=true;
            vector<VariantInfo> &vV = vcfVariants[chrName[ichr]];
            
            std::sort(vV.begin(), vV.end(), 
                [](const VariantInfo &vi1, const VariantInfo &vi2) {
                     return vi1.pos < vi2.pos;
                });
            
                
            //filter: remove variants overlapping deletions
            vector<VariantInfo> vV1;
            vV1.reserve(vV.size());
            uint64 g0=0; //first base after variant
            for (const auto &v : vV) {
                if (v.pos>=g0)
                    vV1.push_back(v);
                g0=max(g0,v.pos+v.seq[0].size());
            };
            P.inOut->logMain << chrName[ichr] <<" filtered out overlapping variants "<< (uint64)vV.size()-(uint64)vV1.size() <<'\n';
            vV=vV1;
            
            for (const auto &v : vV) {
                chrLength1[ichr] += (int64)v.seq[1].size()-(int64)v.seq[0].size();
            };
            P.inOut->logMain << "Transformed chr length difference: " <<chrName[ichr] <<" "<< (int64)chrLength1[ichr]-(int64)chrLength[ichr] <<'\n';
        };
        
        //recalculate chrStart
        chrStart1[0] = 0;
        for (uint32 ichr=0; ichr<chrLength.size(); ichr++) {
            chrStart1[ichr+1]=chrStart1[ichr]+((chrLength1[ichr]+1)/genomeChrBinNbases+1)*genomeChrBinNbases;//+1 makes sure that there is at least one spacer base between chromosomes
            P.inOut->logMain << "Transformed chr start difference: " << chrName[ichr] <<" "<< chrStart1[ichr]-chrStart[ichr] <<'\n';
        };
        
        uint64 nGenome1=chrStart.back();
        P.inOut->logMain << "Old/new genome sizes: " << nGenome <<" "<< nGenome1 <<endl;
        
        //fill Gnew, calculate transformBlocks
        uint64 nG1allocNew;
        char *Gnew=NULL, *G1new=NULL;
        genomeSequenceAllocate(nGenome1, nG1allocNew, Gnew, G1new);
        
        for (uint32 ichr=0; ichr<chrLength.size(); ichr++) {
            
            if (!chtTransformYes[ichr]) {//simple copy
                memcpy(Gnew+chrStart1[ichr], G+chrStart[ichr], chrLength[ichr]);
                transformBlocks.push_back({chrStart[ichr], chrLength[ichr], chrStart1[ichr]});
                continue;
            };
            
            chtTransformYes[ichr]=true;
            vector<VariantInfo> &vV = vcfVariants[chrName[ichr]];
                    
            uint64 iv=0, g1=chrStart1[ichr],  g0=chrStart[ichr];
            transformBlocks.push_back({g0, 0, g1});//first block for the chromosome
            while (g0<chrStart[ichr]+chrLength[ichr]) {
                if (g0==vV[iv].pos-1+chrStart[ichr]) {
                    
                    array<string,2> &seq = vV[iv].seq;
                    
                    //debug
                    char s0[seq[0].size()];
                    convertNucleotidesToNumbers(seq[0].c_str(), s0, seq[0].size());
                    if (memcmp(G+g0, s0, seq[0].size()))
                        cerr <<g0<<" "<<seq[0]<<" "<<G+g0<<endl;
                    //debug
                    
                    char s1[seq[1].size()];
                    convertNucleotidesToNumbers(seq[1].c_str(), s1, seq[1].size());
                    memcpy(Gnew+g1, s1, seq[1].size());
                    g0 += seq[0].size();
                    g1 += seq[1].size();
                    
                    if (vV[iv].len!=0) {//new block
                        //length of the previous block:
                        transformBlocks.back()[1] = g0-seq[0].size() + min(seq[0].size(), seq[1].size()) - transformBlocks.back()[0];
                        transformBlocks.push_back({g0, 0, g1});
                    };
                    
                    if (iv<vV.size()-1) //do not overshoot vV
                        ++iv;
                    
                } else {
                    Gnew[g1]=G[g0];
                    ++g0;
                    ++g1;
                };
            };
            
            if (transformBlocks.back()[1] == 0)
                transformBlocks.back()[1] = g0 - transformBlocks.back()[0];
            
            if (g1!=chrStart1[ichr]+chrLength1[ichr])
                cerr << g1 <<" "<< chrStart1[ichr]+chrLength1[ichr] <<endl;
        };
          
        {//write out transformBlocks
            ofstream & convStream = ofstrOpen(P.pGe.gDir+"/transformGenomeBlocks.tsv",ERROR_OUT, P);
            convStream << transformBlocks.size() <<'\t'<< "-1" <<'\n'; //no -strand
    
            for (auto &tb : transformBlocks) {
                convStream << tb[2] <<'\t'<< tb[1] <<'\t'<< tb[0] <<'\n'; //revert old new for reverse conversion
            };
            convStream.close();
        };
        
        {//transform exonLoci with transformBlocks
            //simple point transformation for start and end
            //start in the gap moves to the right, end in the gap moves to the left
            auto exonLoci1(exonLoci);
            exonLoci1.clear();
            for (auto & exon : exonLoci) {
                
                auto exonS = exon[exS];
                auto exonE = exon[exE];
                auto tBit = std::upper_bound(transformBlocks.begin(), transformBlocks.end(), array<uint64,3> {exonS,0,0},
                               [](const array<uint64,3> &t1, const array<uint64,3> &t2)
                               {
                                  return t1[0] < t2[0];
                               });
                               
               --tBit; //tBit is last block start on the left of exonS
               auto tB=*tBit;
               
               if (exonS < tB[0]+tB[1]) {//exonS inside block
                   exon[exS]=tB[2]+exonS-tB[0];
               } else {
                   exon[exS]=tBit[1][2];//exon start shifts to the next block
               };
               
               while ( exonE > (*tBit)[0]+(*tBit)[1] ) //until exonE is not past the end of the block
                   ++tBit;
               
               tB=*tBit;
               if (exonE >= tB[0]) {//exonE inside block
                   exon[exE]=tB[2]+exonE-tB[0];
               } else {
                   exon[exE]=tBit[-1][2]+tBit[-1][2]-1;//exon end shifts to the end of the previous block
               };
            
               if (exon[exS]<=exon[exE])
                   exonLoci1.push_back(exon);
            };

            P.inOut->logMain << "Transform exons: removed " << exonLoci.size()- exonLoci1.size() <<endl;
            exonLoci=exonLoci1;
            
//             sort(exonLoci.begin(), exonLoci.end(),
//                 [](const array<uint64, exL>& e1, const array<uint64, exL>& e2) {
//                  return e1[exS] < e2[exS];
//                 });
        };
        
        delete[] G1;
        G1=G1new;
        G=Gnew;
        nG1alloc=nG1allocNew;
        nGenome=nGenome1;
        chrStart=chrStart1;
        chrLength=chrLength1;
    };
    
    //diploid: duplicate chromosomes, insert genotypes into sequence
};