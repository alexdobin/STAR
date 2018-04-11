#include "ChimericDetection.h"

ChimericDetection::ChimericDetection(Parameters &Pin, Transcript ***trAll, uint *nWinTr, char** Read1in, Genome &mapGenIn, fstream *ostreamChimJunctionIn, ReadAlign *RAin)
                  : P(Pin), RA(RAin), trAll(trAll), nWinTr(nWinTr), Read1(Read1in), outGen(mapGenIn), ostreamChimJunction(ostreamChimJunctionIn) 
                                    {    
};