#ifndef H_SoloFeatureTypes
#define H_SoloFeatureTypes

namespace SoloFeatureTypes
{
//     enum {Gene=0, GeneFull=1, SJ=2, Transcript3p=3, VelocytoSpliced=4, VelocytoUnspliced=5, VelocytoAmbiguous=6, N=7};
//     const static vector<string> Names={"Gene","GeneFull","SJ","Transcript3p","VelocytoSpliced","VelocytoUnspliced","VelocytoAmbiguous",};
    enum {Gene=0, GeneFull=1, SJ=2, Transcript3p=3, Velocyto=4, N=5};
    const static vector<string> Names={"Gene","GeneFull","SJ","Transcript3p","Velocyto"};
    
};

#endif