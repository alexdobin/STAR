#include "streamFuns.h"
#include "ErrorWarning.h"
#include <fstream>
#ifdef _WIN32
#include <sys/stat.h>
#include <Windows.h>
#else
#include <sys/statvfs.h>
#endif
#include <stdio.h>
#define fstream_Chunk_Max 2147483647

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N) {
    unsigned long long C=0;
    for (unsigned long long ii=0; ii<N/fstream_Chunk_Max; ii++) {
        S.read(A+C,fstream_Chunk_Max);
        C+=S.gcount();
        if (!S.good()) break;
    };
    S.read(A+C,N%fstream_Chunk_Max);
    C+=S.gcount();
    return C;
};

#ifdef _WIN32
int checkDiskSpace(PULARGE_INTEGER pulAvailable, PULARGE_INTEGER pulTotal, PULARGE_INTEGER pulFree){

	if (!GetDiskFreeSpaceEx(NULL, pulAvailable, pulTotal, pulFree))
	{
		DWORD	dwSectorsPerCluster = 0, dwBytesPerSector = 0, dwFreeBytes = 0, dwTotalBytes = 0;
		DWORD	dwFreeClusters = 0, dwTotalClusters = 0;

		if (!GetDiskFreeSpaceA(NULL, &dwSectorsPerCluster,
			&dwBytesPerSector, &dwFreeClusters, &dwTotalClusters)) 
		{
			cout << "Couldn't get free disk space.." << endl;
			return 1;
		}

		dwFreeBytes = dwFreeClusters * dwSectorsPerCluster * dwBytesPerSector;
		dwTotalBytes = dwTotalClusters * dwSectorsPerCluster * dwBytesPerSector;
		// do something with dwFreeBytes & dwTotalBytes
		cout << "Free Bytes " << dwFreeBytes << endl << "Total bytes " << dwTotalBytes << endl;
	}
	return 0; 
}
#endif

void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N, std::string fileName, std::string errorID, Parameters *P) {

#ifndef _WIN32
    struct statvfs statvfsBuf;
    statvfs(fileName.c_str(), &statvfsBuf);
	P->inOut->logMain << "Writing " << N << " bytes into " << fileName << " ; empty space on disk = " << statvfsBuf.f_bavail * statvfsBuf.f_bsize << " bytes ..." << flush;
#else
	ULARGE_INTEGER ulAvailable = { 0 }, ulTotal = { 0 }, ulFree = { 0 };
	checkDiskSpace(&ulAvailable, &ulTotal, &ulFree); 
	P->inOut->logMain << "Writing " << N << " bytes into " << fileName << " ; empty space on disk = " << ulFree.QuadPart << " bytes ..." << flush;
#endif
    
    unsigned long long C=0;
    unsigned long long iC;
    for (iC=0; iC<N/fstream_Chunk_Max; iC++) {
        S.write(A+C,fstream_Chunk_Max);
        C+=fstream_Chunk_Max;
    };
    if (!S.fail()) S.write(A+C,N%fstream_Chunk_Max);
    if (S.fail()) {//failed to write

#ifndef _WIN32
		struct statvfs statvfsBuf;
		statvfs(fileName.c_str(), &statvfsBuf);
#else
		ULARGE_INTEGER ulAvailable = { 0 }, ulTotal = { 0 }, ulFree = { 0 };
		checkDiskSpace(&ulAvailable, &ulTotal, &ulFree);
#endif
		    
//         system(( "ls -lL "+ P->genomeDir + " > "+ P->genomeDir +"/error.info 2>&1").c_str());
//         ifstream error_info((P->genomeDir +"/error.info").c_str());
//         P->inOut->logMain <<error_info.rdbuf();
        
        struct stat statBuf;
        stat(fileName.c_str(), &statBuf);
        
        remove(fileName.c_str());
        
        ostringstream errOut;
        errOut << errorID<<": exiting because of *OUTPUT FILE* error: could not write the output file "<< fileName <<"\n";
        errOut << "fail()=" <<S.fail() <<" ; bad()="<< S.bad()<<"\n";
        errOut << "Error while trying to write chunk # " << iC << "; "<< C << " bytes\n";
        errOut << "File size full = "<< N <<" bytes\n";        
        errOut << "File size on disk = " << statBuf.st_size<<" bytes\n";
        errOut << "Solution: check that you have enough space on the disk\n";
#ifndef _WIN32
        errOut << "Empty space on disk = " << statvfsBuf.f_bavail * statvfsBuf.f_bsize <<" bytes\n";
#else
		errOut << "Empty space on disk = " << ulFree.QuadPart << " bytes\n";
#endif
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_FILE_WRITE, *P);
    };
    P->inOut->logMain << " done\n" <<flush;
};

std::ofstream & ofstrOpen (std::string fileName, std::string errorID, Parameters *P) {//open file 'fileName', generate error if cannot open
    std::ofstream & ofStream = *new std::ofstream(fileName.c_str(), std::fstream::out | std::fstream::trunc | std::fstream::binary);
    if (ofStream.fail()) {//
//         dir1=fileName.substr(0,fileName.find_last_of("/")+1);
//         if (dir1=="") dir1="./";
        ostringstream errOut;
        errOut << errorID<<": exiting because of *OUTPUT FILE* error: could not create output file "<< fileName <<"\n";
        errOut << "Solution: check that the path exists and you have write permission for this file\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_FILE_OPEN, *P);
    };    
    return ofStream;
};


std::ifstream & ifstrOpen (std::string fileName, std::string errorID, std::string solutionString, Parameters *P) {
    //open file 'fileName', generate error if cannot open
    std::ifstream & ifStream = *new std::ifstream(fileName.c_str(), ios_base::in | ios_base::binary);
    if (ifStream.fail()) {//
        ostringstream errOut;
        errOut << errorID<<": exiting because of *INPUT FILE* error: could not open input file "<< fileName <<"\n";
        errOut << "Solution: check that the file exists and you have read permission for this file\n";
        if (solutionString.size()>0) {
            errOut << "          "<< solutionString <<"\n";
        };
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_FILE_OPEN, *P);
    };    
    return ifStream;
};

ifstream & ifstrOpenGenomeFile (std::string fileName, std::string errorID, Parameters *P) {
     //open one of the genome files
     return ifstrOpen(P->genomeDir+"/"+fileName, errorID,  "if this file is missing from the genome directory, you will need to *re-generate the genome*", P);
};

void copyFile(string fileIn, string fileOut)
{//copy fileIn into FileOut
    std::ifstream  src(fileIn, std::ios::binary);
    std::ofstream  dst(fileOut,   std::ios::binary);
    dst << src.rdbuf();
};
