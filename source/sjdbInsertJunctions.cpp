#include "sjdbInsertJunctions.h"
#include "sjdbLoadFromStream.h"
#include "sjdbLoadFromFiles.h"
#include "sjdbPrepare.h"
#include "ErrorWarning.h"
#include "loadGTF.h"
#include "sjdbBuildIndex.h"
#include "streamFuns.h"
#include "genomeParametersWrite.h"

void sjdbInsertJunctions(Parameters * P, Parameters * P1, Genome & genome, SjdbClass & sjdbLoci) 
{
    time_t rawtime;

    if (P->sjdbN>0 && sjdbLoci.chr.size()==0)
    {//load from the saved genome, only if the loading did not happen already (if sjdb insertion happens at the 1st pass, sjdbLoci will be populated
        ifstream & sjdbStreamIn = ifstrOpen(P->genomeDir+"/sjdbList.out.tab", ERROR_OUT, "SOLUTION: re-generate the genome in genomeDir=" + P->genomeDir, P);
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        sjdbLoci.priority.resize(sjdbLoci.chr.size(),30);
        time ( &rawtime );
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the generated genome " << P->genomeDir+"/sjdbList.out.tab" <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    };    
    
    if (P->twoPass.pass2)
    {//load 1st pass new junctions
     //sjdbLoci already contains the junctions from before 1st pass
        ifstream sjdbStreamIn ( P->twoPass.pass1sjFile.c_str() );   
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL INPUT error, could not open input file with junctions from the 1st pass=" << P->twoPass.pass1sjFile <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
        sjdbLoadFromStream(sjdbStreamIn, sjdbLoci);
        sjdbLoci.priority.resize(sjdbLoci.chr.size(),0);
        time ( &rawtime );
        P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the 1st pass file: " << P->twoPass.pass1sjFile <<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
    } else 
    {//loading junctions from GTF or tab or from the saved genome is only allowed at the 1st pass
     //at the 2nd pass these are already in the sjdbLoci
        
        if (P->sjdbFileChrStartEnd.at(0)!="-")
        {//load from junction files
            sjdbLoadFromFiles(P,sjdbLoci);
            sjdbLoci.priority.resize(sjdbLoci.chr.size(),10);
            time ( &rawtime );        
            P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the sjdbFileChrStartEnd file(s), " << sjdbLoci.chr.size()<<" total junctions\n\n";        
        };

        if (P->sjdbGTFfile!="-")
        {//load from GTF
            loadGTF(sjdbLoci, P, P->sjdbInsert.outDir);
            sjdbLoci.priority.resize(sjdbLoci.chr.size(),20);
            time ( &rawtime );
            P->inOut->logMain << timeMonthDayTime(rawtime) << "   Loaded database junctions from the GTF file: " << P->sjdbGTFfile<<": "<<sjdbLoci.chr.size()<<" total junctions\n\n";
        };
    };
    
    char *Gsj=new char [2*P->sjdbLength*sjdbLoci.chr.size()+1];//array to store junction sequences, will be filled in sjdbPrepare
    sjdbPrepare (sjdbLoci, P, P->chrStart[P->nChrReal], P->sjdbInsert.outDir, genome.G, Gsj);//P->nGenome - change when replacing junctions
    time ( &rawtime );
    P->inOut->logMain  << timeMonthDayTime(rawtime) << "   Finished preparing junctions" <<endl;

    if (P->sjdbN>P->limitSjdbInsertNsj)
    {
        ostringstream errOut;
        errOut << "Fatal LIMIT error: the number of junctions to be inserted on the fly ="<<P->sjdbN<<" is larger than the limitSjdbInsertNsj="<<P->limitSjdbInsertNsj<<"\n";                errOut << "Fatal LIMIT error: the number of junctions to be inserted on the fly ="<<P->sjdbN<<" is larger than the limitSjdbInsertNsj="<<P->limitSjdbInsertNsj<<"\n";
        errOut << "SOLUTION: re-run with at least --limitSjdbInsertNsj "<<P->sjdbN<<"\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };
    //insert junctions into the genome and SA and SAi
    sjdbBuildIndex (P, P1, Gsj, genome.G, genome.SA, (P->twoPass.pass2 ? genome.SApass2 : genome.SApass1), genome.SAi);
    delete [] Gsj; //junction sequences have been added to G
    time ( &rawtime ); 
    P->inOut->logMain     << timeMonthDayTime(rawtime) << " ..... Finished inserting junctions into genome" <<endl;    

    if (P->sjdbInsert.save=="All")
    {//save and copy all genome files into sjdbInsert.outDir, except those created above
        if (P->genomeDir != P->sjdbInsert.outDir)
        {
            copyFile(P->genomeDir+"/chrName.txt", P->sjdbInsert.outDir+"/chrName.txt");
            copyFile(P->genomeDir+"/chrStart.txt", P->sjdbInsert.outDir+"/chrStart.txt");
            copyFile(P->genomeDir+"/chrNameLength.txt", P->sjdbInsert.outDir+"/chrNameLength.txt");
            copyFile(P->genomeDir+"/chrLength.txt", P->sjdbInsert.outDir+"/chrLength.txt");
        };
        
        genomeParametersWrite(P->sjdbInsert.outDir+("/genomeParameters.txt"), P, ERROR_OUT);
        
        ofstream & genomeOut = ofstrOpen(P->sjdbInsert.outDir+"/Genome",ERROR_OUT, P);
        fstreamWriteBig(genomeOut,genome.G,P->nGenome,P->sjdbInsert.outDir+"/Genome",ERROR_OUT,P);
        genomeOut.close();

        ofstream & saOut = ofstrOpen(P->sjdbInsert.outDir+"/SA",ERROR_OUT, P);
        fstreamWriteBig(saOut,(char*) genome.SA.charArray, (streamsize) genome.SA.lengthByte, P->sjdbInsert.outDir+"/SA",ERROR_OUT,P);
        saOut.close();

        ofstream & saIndexOut = ofstrOpen(P->sjdbInsert.outDir+"/SAindex",ERROR_OUT, P);
        fstreamWriteBig(saIndexOut, (char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases),P->sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        fstreamWriteBig(saIndexOut, (char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1),P->sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        fstreamWriteBig(saIndexOut,  genome.SAi.charArray, genome.SAi.lengthByte,P->sjdbInsert.outDir+"/SAindex",ERROR_OUT,P);
        saIndexOut.close();
    };

    //re-calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;
};
