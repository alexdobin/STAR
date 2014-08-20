#ifndef CODE_signalFromBAM
#define CODE_signalFromBAM
#include "htslib/htslib/sam.h"
#include  <fstream>
#include <string>
#include "Stats.h"

using namespace std;

void signalFromBAM(const string bamFileName, const string sigFileName, const bool flagStranded, const int signalType, const string outWigReferences);

#endif
