#ifndef CODE_GenomeInsertAnnotations
#define CODE_GenomeInsertAnnotations

#include "IncludeDefine.h"

class Parameters;

void genomeInsertCopyAnnotationSidecars(const string &dirIn, const string &dirOut);
void genomeInsertMergeAnnotationSidecars(const string &dirBase, const string &dirInsert, const string &dirOut, Parameters &P, bool copySjdbFiles);

#endif
