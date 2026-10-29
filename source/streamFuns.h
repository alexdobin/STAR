#ifndef CODE_streamFuns
#define CODE_streamFuns

#include "Parameters.h"
#include <fstream>

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N);
unsigned long long fstreamReadBigParallel(std::ifstream &S, const std::string &fileName, unsigned long long offset, char* A, unsigned long long N, int nThreads);//pread N bytes from offset with nThreads threads; falls back to fstreamReadBig(S,...), which needs S at offset
void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N, std::string fileName, std::string errorID, Parameters &P) ;

fstream  &fstrOpen  (std::string fileName, std::string errorID, Parameters &P, bool flagDelete);
ofstream &ofstrOpen (std::string fileName, std::string errorID, Parameters &P);
ifstream &ifstrOpen (std::string fileName, std::string errorID, std::string solutionString, Parameters &P);
ifstream &ifstrOpenGenomeFile (std::string fileName, std::string errorID, Parameters &P);
void createDirectory(const string dirPathIn, const mode_t dirPerm, const string dirParameter, Parameters &P);

void copyFile(string fileIn, string fileOut);
#endif
