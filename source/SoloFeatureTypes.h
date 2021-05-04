#ifndef H_SoloFeatureTypes
#define H_SoloFeatureTypes

namespace SoloFeatureTypes
{
//     enum {Gene=0, GeneFull=1, SJ=2, Transcript3p=3, VelocytoSpliced=4, VelocytoUnspliced=5, VelocytoAmbiguous=6, N=7};
//     const static vector<string> Names={"Gene","GeneFull","SJ","Transcript3p","VelocytoSpliced","VelocytoUnspliced","VelocytoAmbiguous",};
    enum {SJ=0, Transcript3p=1, GeneFull=2, Gene=3, VelocytoSimple=4, Velocyto=5, N=6};
    const static vector<string> Names={"SJ", "Transcript3p", "GeneFull", "Gene", "VelocytoSimple", "Velocyto"};
};

#endif