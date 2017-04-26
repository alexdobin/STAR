#ifndef CODE_streamFuns
#define CODE_streamFuns

#include "Parameters.h"
#include <fstream>

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N);
void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N, std::string fileName, std::string errorID, Parameters &P) ;

ofstream & ofstrOpen (std::string fileName, std::string errorID, Parameters &P);
ifstream & ifstrOpen (std::string fileName, std::string errorID, std::string solutionString, Parameters &P);
ifstream & ifstrOpenGenomeFile (std::string fileName, std::string errorID, Parameters &P);

void copyFile(string fileIn, string fileOut);
#endif
