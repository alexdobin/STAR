#ifndef CODE_ParametersChimeric
#define CODE_ParametersChimeric

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"

class ParametersChimeric
{//
    public:
        uint chimSegmentMin, chimJunctionOverhangMin; //min chimeric donor/acceptor length
        uint chimSegmentReadGapMax; //max read gap for stitching chimeric windows
        int chimScoreMin,chimScoreDropMax,chimScoreSeparation, chimScoreJunctionNonGTAG; //min chimeric score
        uint chimMainSegmentMultNmax;
        vector <string> chimFilter;

        struct
        {
            bool genomicN;
        } filter;
        struct
        {
            vector <string> type;
            bool bam;
            bool bamHardClip;
        } out;
    private:
};

#endif