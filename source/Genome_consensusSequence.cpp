#include "Genome.h"
#include "SequenceFuns.h"
#include "streamFuns.h"

void Genome::consensusSequence() {
    if (pGe.gConsensusFile!="-") {//load consensus SNPs
        ifstream &consIn=ifstrOpen(pGe.gConsensusFile, ERROR_OUT, "SOLUTION: check path and permission for the --genomeConsensusFile file" + pGe.gConsensusFile, P);

        map<string,uint> chrStartMap;
        for (uint ii=0;ii<nChrReal;ii++) {
            chrStartMap.insert(std::pair <string,uint> (chrName[ii], chrStart[ii]));
        };

        uint nInserted=0, nWrongChr=0, nWrongRef=0, nRefN=0;
        while (consIn.good()) {
            string chr1, refIn, altIn, dummy;
            uint start1;
            char ref1,alt1;

            consIn >> chr1 >> start1 >> dummy >> refIn >> altIn;
            consIn.ignore(numeric_limits<streamsize>::max(),'\n');

            convertNucleotidesToNumbers(refIn.c_str(),&ref1,1);
            convertNucleotidesToNumbers(altIn.c_str(),&alt1,1);
            --start1;//VCF positions are 1-based

            if (chrStartMap.count(chr1)==1) {//otherwise just skip
                start1+=chrStartMap[chr1];
                if (G[start1]>3)
                    ++nRefN;

                if (G[start1]==ref1 || G[start1]>3) {
                    G[start1]=alt1;
                    ++nInserted;
                } else {
                    ++nWrongRef;
                    P.inOut->logMain << "WARNING: reference allele in consensus file does not agree with reference genome base: ";
                    P.inOut->logMain << chr1 <<"   "<< start1-chrStartMap[chr1] <<"   "<< (int) G[start1]<<"   "<< (int) ref1<<"   "<< (int) alt1<<"\n";
                };
            } else {
                ++nWrongChr;
            };
        };
        P.inOut->logMain <<"Inserted consensus variants: " << nInserted <<", including reference N-base:"<< nRefN <<", wrong chromosome: " << nWrongChr<< ", wrong reference base: " << nWrongRef << endl;
    };
};