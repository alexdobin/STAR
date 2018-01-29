#ifndef CODE_ParametersChimeric
#define CODE_ParametersChimeric

#include "IncludeDefine.h"

class ParametersChimeric
{//
    public:
        uint segmentMin, junctionOverhangMin; //min chimeric donor/acceptor length
        uint segmentReadGapMax; //max read gap for stitching chimeric windows
        int scoreMin,scoreDropMax,scoreSeparation, scoreJunctionNonGTAG; //min chimeric score
        uint mainSegmentMultNmax;
        
        uint multimapScoreRange, multimapNmax, nonchimScoreDropMin;

        struct
        {
            vector <string> stringIn;
            bool genomicN;
        } filter;

        struct
        {
            vector <string> type;
            bool bam;
            bool bamHardClip;
            bool samOld;
            bool junctions;
        } out;
    private:
};

#endif