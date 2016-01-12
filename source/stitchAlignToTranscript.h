#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

#ifdef __cplusplus
extern "C" {
#endif

intScore stitchAlignToTranscript(uint rAend, uint gAend, uint rBstart, uint gBstart, uint L, uint iFragB, uint sjAB, Parameters* P, char* R, char* Q, char* G,  Transcript *trA, uint outFilterMismatchNmaxTotal);

#ifdef __cplusplus
}
#endif