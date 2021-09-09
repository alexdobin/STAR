#ifndef H_SoloFeatureTypes
#define H_SoloFeatureTypes

namespace SoloFeatureTypes
{
//     enum {Gene=0, GeneFull=1, SJ=2, Transcript3p=3, VelocytoSpliced=4, VelocytoUnspliced=5, VelocytoAmbiguous=6, N=7};
//     const static vector<string> Names={"Gene","GeneFull","SJ","Transcript3p","VelocytoSpliced","VelocytoUnspliced","VelocytoAmbiguous",};
    enum {SJ=0, Transcript3p=1, GeneFull=2, GeneFull_ExonOverIntron=3, GeneFull_Ex50pAS=4, Gene=5, VelocytoSimple=6, Velocyto=7, N=8};
    const static vector<string> Names={"SJ", "Transcript3p", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS", "Gene", "VelocytoSimple", "Velocyto"};
};

#endif