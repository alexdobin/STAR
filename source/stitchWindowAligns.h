#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "extendAlign.h"
#include "stitchAlignToTranscript.h"
#include "ReadAlign.h"

void stitchWindowAligns(uint iA, uint nA, int Score, bool WAincl[], uint tR2, uint tG2, Transcript trA, \
                        uint Lread, uiWA* WA, char* R, Genome &mapGen, \
                        Parameters& P, Transcript** wTr, uint* nWinTr, ReadAlign *RA);
    //recursively stitch aligns for one gene
    //*nWinTr - number of transcripts for the current window
