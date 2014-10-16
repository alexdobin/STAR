#ifndef CODE_streamFuns
#define CODE_streamFuns

#include "Parameters.h"
#include <fstream>

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N);
void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N);

void openOfstream (std::string fileName, std::string errorID, Parameters *P, ofstream & ofStream);

#endif
