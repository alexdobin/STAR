#include "readLoad.h"
#include "ErrorWarning.h"

static inline bool isSpaceClassic(char c)
{//isspace in the classic locale, which istream >> uses here
    return c==' ' || (c>='\t' && c<='\r');
};

static bool readNameLineParse(const string &line1, char* readName, uint &iReadAll, char &readFilter, uint32 &readFilesIndex, string &readNameExtra)
{//the line processChunks writes, "name iReadAll readFilter readFilesIndex[ rest]", parsed as readLoad's istringstream would;
 //false, changing nothing, for any other line (a sign, a long number or name, a missing field)
    const char *p=line1.data(), *pEnd=p+line1.size();

    while (p<pEnd && isSpaceClassic(*p))
        ++p;
    const char *name1=p;
    while (p<pEnd && !isSpaceClassic(*p))
        ++p;
    uint nameLength=p-name1;
    if (nameLength==0 || nameLength>=DEF_readNameLengthMax-1)
        return false;

    while (p<pEnd && isSpaceClassic(*p))
        ++p;
    const char *num1=p;
    uint iReadAll1=0;
    while (p<pEnd && *p>='0' && *p<='9' && p-num1<19) {//19 digits always fit in uint
        iReadAll1=iReadAll1*10+(*p-'0');
        ++p;
    };
    if (p==num1 || p==pEnd || (*p>='0' && *p<='9'))
        return false;

    while (p<pEnd && isSpaceClassic(*p)) //>> char skips whitespace, then takes any one char
        ++p;
    if (p==pEnd)
        return false;
    char readFilter1=*p;
    ++p;

    while (p<pEnd && isSpaceClassic(*p))
        ++p;
    const char *num2=p;
    uint32 readFilesIndex1=0;
    while (p<pEnd && *p>='0' && *p<='9' && p-num2<9) {//9 digits always fit in uint32
        readFilesIndex1=readFilesIndex1*10+(*p-'0');
        ++p;
    };
    if (p==num2 || (p<pEnd && !isSpaceClassic(*p)))
        return false;

    memcpy(readName, name1, nameLength);
    readName[nameLength]=0;
    iReadAll=iReadAll1;
    readFilter=readFilter1;
    readFilesIndex=readFilesIndex1;

    while (p<pEnd && isSpaceClassic(*p))
        ++p;
    if (p<pEnd) //at the end of the line, >> std::ws and getline fail and leave readNameExtra as it was
        readNameExtra.assign(p, pEnd-p);
    return true;
};

int readLoad(istream& readInStream, Parameters& P, uint& Lread, uint& LreadOriginal, \
            char* readName, char* Seq, char* SeqNum, char* Qual, vector<ClipMate> &clipOneMate, \
            uint &iReadAll, uint32 &readFilesIndex, char &readFilter, string &readNameExtra)
{//load one read from a stream
    int readFileType=0;

    if (readInStream.peek()!='@' && readInStream.peek()!='>') 
        return -1; //end of the stream

    string line1;
    
    getline(readInStream, line1); //read 1st line
    if (!readNameLineParse(line1, readName, iReadAll, readFilter, readFilesIndex, readNameExtra)) {
        istringstream line1str(line1);

        readName[0]=0;//clear char array
        line1str >> readName; //TODO check that it does not overflow the array
        if (strlen(readName)>=DEF_readNameLengthMax-1) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: read name is too long:" << readInStream.gcount()<<"\nRead Name="<<readName<<"\nDEF_readNameLengthMax="<<DEF_readNameLengthMax<<'\n';
            errOut << "SOLUTION: increase DEF_readNameLengthMax in IncludeDefine.h and re-compile STAR\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };  

        line1str >> iReadAll >> readFilter >> readFilesIndex; //extract read number TODO this is needed only for one mate
        line1str >> std::ws; //skip whitespace before readNameExtra
        getline(line1str, readNameExtra);
    };

    readInStream.getline(Seq,DEF_readSeqLengthMax+1); //extract sequence

    Lread=readInStream.gcount()-1;

    if (Lread<1) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: short read sequence line: " << Lread <<"\n";
        errOut << "Read Name="<< readName <<'\n'<< "Read Sequence=\"" << Seq <<"\"\nDEF_readNameLengthMax="<< DEF_readNameLengthMax <<'\n'<< "DEF_readSeqLengthMax="<<DEF_readSeqLengthMax<<'\n';
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    if (Lread>DEF_readSeqLengthMax) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: Lread>=" << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<'\n'<< "Read Name="<<readName<<'\n';
        errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    LreadOriginal=Lread;

    convertNucleotidesToNumbers(Seq,SeqNum,Lread);
    
    if (readName[0]=='@') {//fastq format: read the '+' line
        readFileType=2;
        readInStream.get(clipOneMate[0].clippedInfo); //this is used if clipChunk is activated, only 5' for now
        readInStream.ignore(); //ignore one char: \n
    };

    //cout << readName <<'\t' << (uint32) line3char << '\n';
    
    clipOneMate[0].clip(Lread, SeqNum); //5p clip
    clipOneMate[1].clip(Lread, SeqNum); //3p clip

    if (readName[0]=='@') {//fastq format, read qualities
        readInStream.getline(Qual,DEF_readSeqLengthMax);//read qualities
        if ((uint) readInStream.gcount() != LreadOriginal+1) {//inconsistent read sequence and quality
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length\n";
            errOut << readName <<'\n'<< Seq <<'\n'<< Qual <<'\n'<< "SOLUTION: fix your fastq file\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        if (P.outQSconversionAdd!=0) {
            for (uint ii=0;ii<LreadOriginal;ii++) {
                int qs=int(Qual[ii])+P.outQSconversionAdd;
                if (qs<33) {
                    qs=33;
                } else if (qs>126) {
                    qs=126;
                };
                Qual[ii]=qs;
            };
        };
    
    } else if (readName[0]=='>') {//fasta format, assign Qtop to all qualities
        readFileType=1;
        for (uint ii=0;ii<LreadOriginal;ii++) 
            Qual[ii]='A';
        Qual[LreadOriginal]=0;
    } else {//wrong header
        ostringstream errOut;
        errOut <<"Unknown reads file format: header line does not start with @ or > : "<< readName<<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    //trim read name TODO this is needed only for one mate
    for (uint ii=0; ii<P.readNameSeparatorChar.size(); ii++) {
        char* pSlash=strchr(readName,P.readNameSeparatorChar.at(ii)); //trim everything after ' '
        if (pSlash!=NULL) *pSlash=0;
    };
    return readFileType;
};
