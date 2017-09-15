#include "Chain.h"
#include "streamFuns.h"
#include "serviceFuns.cpp"

Chain::Chain(Parameters &Pin, string chainFileNameIn) : P(Pin), chainFileName(chainFileNameIn)
{
    chainLoad();
};

void Chain::chainLoad()
{
    ifstream &streamIn = ifstrOpen(chainFileName, ERROR_OUT, "SOLUTION: check path and permission for the chain file" + chainFileName, P);
      
    string chr1;//current chromsome 1 (old)
    
    while (streamIn.good())
    {        
        string line1;
        getline(streamIn,line1);
        istringstream line1str(line1);

        vector <string>  fields(13);

        for (int ii=0;ii<4;ii++)
            line1str >> fields[ii];
        if (fields[0]=="")
        {//empty line, continue
        } else if (fields[1]=="")
        {//end of chain
            chrChains[chr1].bLen.push_back(std::stoi(fields[0]));//read the last block length
            chrChains[chr1].bN=chrChains[chr1].bLen.size();
        } else if (fields[3]=="")
        {//normal chain block
            chrChains[chr1].bLen.push_back(std::stoi(fields[0]));
            
            uint s=chrChains[chr1].bStart1.back() + chrChains[chr1].bLen.back() + std::stoi(fields[1]);//prev start + length + shift
            chrChains[chr1].bStart1.push_back(s);
            
            s=chrChains[chr1].bStart2.back() + chrChains[chr1].bLen.back() + std::stoi(fields[2]);//prev start + length + shift
            chrChains[chr1].bStart2.push_back(s);
        } else
        {//chain header
            //chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
            //    0     1     2     3     4       5      6    7     8     9       10   11 12 
            
            for (int ii=4;ii<13;ii++)
                line1str >> fields[ii]; //read all the fields
            
            chr1=fields[2];
            chrChains[chr1].chr1=chr1;
            chrChains[chr1].chr2=fields[7];//NOTE: the whole procedure (for now) only works for single chain per chromosome
            chrChains[chr1].bStart1.push_back(std::stoi(fields[5]));
            chrChains[chr1].bStart2.push_back(std::stoi(fields[10]));
        };
    };
};

void Chain::liftOverGTF(string gtfFileName, string outFileName)
{//simple arithmetic lift-over of the GTF file
    ifstream &streamIn = ifstrOpen(gtfFileName, ERROR_OUT, "SOLUTION: check path and permission for the GTF file" + gtfFileName, P);
    ofstream &streamOut = ofstrOpen(outFileName, ERROR_OUT, P);
    ofstream &streamOutUnlifted = ofstrOpen(outFileName+".unlifted", ERROR_OUT, P);
    
    while (streamIn.good())
    {
        string line1;
        getline(streamIn,line1);
        istringstream line1str(line1);

        string chr1;
        line1str >> chr1;
        
        if (chr1=="" || chr1.substr(0,1)=="#")
            continue;//empty or comment line
        
        if (chrChains.count(chr1)==0)
            exitWithError("GTF contains chromosome " + chr1 + " not present in the chain file " + chainFileName,std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        
        OneChain *ch1 = & chrChains[chr1];//the chain for the chr1
               
        string str1,str2;
        line1str >> str1 >> str2;//fields 2,3
        
        uint c1, c2[2]; //coordinates: 1/2 (old/new)
        
        for (int ii=0;ii<2;ii++)
        {//read and transform the two coordinates
            line1str >> c1;
            int32 i1 = binarySearch1a <uint> (c1, ch1->bStart1.data(), ch1->bN);
            
            c2[ii]=-1;//-1 means impossible to lift this end
           
            if (i1>=0 && c1 < ch1->bStart1[i1]+ch1->bLen[i1])
            {//c1 is inside the block, simple transformation
                c2[ii]=ch1->bStart2[i1] + c1 - ch1->bStart1[i1];
            } else
            {//c1 is outside of the block 
                if (ii==0 && i1 < (int32) ch1->bN-1)
                {//left end => c2 will be at the start of the next block
                        c2[ii]=ch1->bStart2[i1+1]; //if i1=-1, it will work = start of the 0-tn blocl
                } else if (ii==1 && i1 >= 0)
                {
                    c2[ii]=ch1->bStart2[i1]+ch1->bLen[i1]-1;
                };
            };
        };
        if (c2[0]!=-1llu && c2[1]!=-1llu && c2[1]>=c2[0])
        {//good conversion
            streamOut << ch1->chr2 <<"\t"<< str1 <<"\t"<< str2 <<"\t"<<c2[0] <<"\t"<< c2[1] << line1str.rdbuf() << "\n";
        } else {
            streamOutUnlifted << line1 <<"\n";
        };
    };
    streamOutUnlifted.close();
    streamOut.close();
};

//             if (mapGen.chrNameIndex.count(oldname))
//             {
//                 ostringstream errOut;
//                 errOut << "EXITING because of fatal INPUT file error: chain file contains chromosome (scaffold) not present in the genome " << oldname <<"\n";
//                 errOut << ERROR_OUT << "\n";
//                 exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
//             };            
//             uint ichr=mapGen.chrNameIndex[oldname];//chr index in the genome list
//             bStart1[bN] += mapGen.chrLength[ichr];//whole genome chain - shift by chr start
