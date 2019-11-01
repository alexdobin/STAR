#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Transcript.h"
#include "extendAlign.h"
#include "stitchAlignToTranscript.h"

void stitchWindowAligns(uint iA, uint nA, int Score, bool WAincl[], uint tR2,
                        uint tG2, const Transcript &trA, uint Lread, const uiWA *WA,
                        const char *R, const Genome &mapGen, const Parameters &P,
                        Transcript **wTr, uint *nWinTr, ReadAlign *RA);
// recursively stitch aligns for one gene
//*nWinTr - number of transcripts for the current window
