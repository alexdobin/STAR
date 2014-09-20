#ifndef CODE_signalFromBAM
#define CODE_signalFromBAM
#include "htslib/htslib/sam.h"
#include  <fstream>
#include <string>
#include "Stats.h"
#include "Parameters.h"

using namespace std;

void signalFromBAM(const string bamFileName, const string sigFileName, Parameters P);

#endif
