#ifndef INCLUDEDEFINE_DEF
#define INCLUDEDEFINE_DEF

//standard libs
#include <algorithm>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctime>
#include <iomanip>
#include <vector>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <limits>
#include <stdint.h>
#include <omp.h>

#include "VERSION"

#define ERROR_OUT string ( __FILE__ ) +":"+ to_string ( (uint) __LINE__ ) +":"+ string ( __FUNCTION__ )

//external libs
#define SAMTOOLS_BGZF_H "htslib/htslib/bgzf.h"
#define SAMTOOLS_SAM_H  "htslib/htslib/sam.h"

using namespace std;

#ifdef COMPILE_FOR_MAC
  //some Mac's idiosyncrasies: standard SHM libraries are very old and missing some definitions
  #define SHM_NORESERVE 0
#endif

#if defined(__mips__) && !defined(SHM_NORESERVE)
#define SHM_NORESERVE 010000
#endif

typedef int8_t int8;
typedef uint8_t uint8;

#define uint unsigned long long
#define sint signed long long
#define uint64 unsigned long long
#define uint32 unsigned int
#define uint16 unsigned short int
#define uchar unsigned char
#define int64 long long
#define int32 int

// this is gcc extension, may need to redefine for other compilers
#define uint128 __uint128_t

#define GENOME_spacingChar 5

#define uintWinBin unsigned short
#define uintWinBinMax numeric_limits<uint16>::max()


#define intSWscore int
#define intScore int

#define scoreMatch 1


//cleaned
//output
#define BAMoutput_oneAlignMaxBytes 100000


//SAM attributes
#define ATTR_NH 1
#define ATTR_HI 2
#define ATTR_AS 3
#define ATTR_NM 4
#define ATTR_MD 5
#define ATTR_nM 6
#define ATTR_jM 7
#define ATTR_jI 8
#define ATTR_XS 9
#define ATTR_RG 10
#define ATTR_vG 11
#define ATTR_vA 12
#define ATTR_vW 13
#define ATTR_ch 14
#define ATTR_MC 15
#define ATTR_rB 16


//BAM definitions
#define BAM_CIGAR_MaxSize 10000
#define BAM_CIGAR_OperationShift 4
#define BAM_CIGAR_M 0
#define BAM_CIGAR_I 1
#define BAM_CIGAR_D 2
#define BAM_CIGAR_N 3
#define BAM_CIGAR_S 4
#define BAM_CIGAR_H 5
#define BAM_CIGAR_P 6
#define BAM_CIGAR_EQ 7
#define BAM_CIGAR_X 8


#if defined COMPILE_FOR_LONG_READS
    #define MAX_N_EXONS 1000
    #define BAM_ATTR_MaxSize 10000
#else
    #define MAX_N_EXONS 20
    #define BAM_ATTR_MaxSize 1000
#endif

//input reads
#define MAX_N_MATES 2
#define DEF_readNameLengthMax 50000
#if defined COMPILE_FOR_LONG_READS
    #define DEF_readSeqLengthMax 500000
#else
    #define DEF_readSeqLengthMax 650
#endif

#if (DEF_readNameLengthMax > DEF_readSeqLengthMax)
        #define DEF_readNameSeqLengthMax DEF_readNameLengthMax
#else
        #define DEF_readNameSeqLengthMax DEF_readSeqLengthMax
#endif

#define EXIT_CODE_BUG 101
#define EXIT_CODE_PARAMETER 102
#define EXIT_CODE_RUNTIME 103
#define EXIT_CODE_INPUT_FILES 104
#define EXIT_CODE_GENOME_FILES 105
#define EXIT_CODE_SHM 106
#define EXIT_CODE_GENOME_LOADING_WAITED_TOO_LONG 107
#define EXIT_CODE_MEMORY_ALLOCATION 108
#define EXIT_CODE_FILE_OPEN 109
#define EXIT_CODE_FILE_WRITE 110
#define EXIT_CODE_INCONSISTENT_DATA 111

//cleaned-end


//exit codes
#define EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS 101

#define SJ_MOTIF_SIZE 7 //number of recorded SJ motifs
#define SJ_SAM_AnnotatedMotifShift 20

#define EXTEND_ORDER 1 //1-first extend to the 5' of the read, then 3'; 2- first extend to the left, then to the right

#define MAX_N_FRAG 2
#define MARK_FRAG_SPACER_BASE 11
#define MAX_N_CHIMERAS 5
#define MAX_N_MULTMAP 100000 //max number of multiple mappers
#define MAX_SJ_REPEAT_SEARCH 255 //max length of a repeat to search around a SJ
#define MAX_QS_VALUE 60
#define MAX_OUTPUT_FLAG 10

#define PC_rStart 0
#define PC_Length 1
#define PC_Str 2
#define PC_Dir 3
#define PC_Nrep 4
#define PC_SAstart 5
#define PC_SAend 6
#define PC_iFrag 7
#define PC_SIZE 8

#define WC_Str 0
#define WC_Chr 1
#define WC_gStart 2
#define WC_gEnd 3
#define WC_SIZE 4

#define WA_Length 0
#define WA_rStart 1
#define WA_gStart 2
#define WA_Nrep 3
#define WA_Anchor 4
#define WA_iFrag 5
#define WA_sjA 6
#define WA_SIZE 7

#define EX_R 0
#define EX_G 1
#define EX_L 2
#define EX_iFrag 3
#define EX_sjA 4
#define EX_SIZE 5

//mapType
#define MT_PE 0 //paired end type
#define MT_SIZE 5

#define MARKER_ALL_PIECES_EXCEED_seedMultimapNmax 999901 //marks the reads that map too many time, more than seedMultimapNmax
#define MARKER_NO_UNIQUE_PIECES 999902 //the best transcript does not contain any unique pieces
#define MARKER_NO_GOOD_WINDOW 999903 //did not find any good windows
#define MARKER_NO_GOOD_PIECES 999904
#define MARKER_TOO_MANY_ANCHORS_PER_WINDOW 999905
#define MARKER_MAX_N_MULT_EXCEEDED 999906
#define MARKER_FULL_LENGTH_MULTIMAPPER_EXCEEDED_alignWindowsPerReadNmax 999907
#define MARKER_ALL_PIECES_EXCEEDED_winAnchorMultimapNmax 999908
#define MARKER_TOO_MANY_CHIMERAS 999909
#define MARKER_READ_TOO_SHORT 999910

#define PEMARKER_SINGLE_END 0
#define PEMARKER_PAIR 1
#define PEMARKER_ONE_END 3
#define PEMARKER_TOO_MANY_PAIRS 5
#define PEMARKER_CHIMERIC_PAIRS 7
#define PEMARKER_CHIMERIC_SJ_READ1 221
#define PEMARKER_CHIMERIC_SJ_READ2 223
#define PEMARKER_CHIMERIC_SJ_READ1and2 225
#define PEMARKER_SINGLE_END_NOTMAPPED 1001


typedef uint uiPC[PC_SIZE];
typedef uint uiWC[WC_SIZE];
typedef uint uiWA[WA_SIZE];

// debugging
//#define DEBUG_Nread 1000000
//#define DEBUG
#if defined DEBUG
    #define DEBUG_stitch
    #define DEBUG_Nread 200000
    #define DEBUG_NreadStart 1
    #define DEBUG_extend
#endif

// #define DEBUG_NreadStart 500000

#endif
