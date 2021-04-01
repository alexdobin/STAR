#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdio>
#include <limits>
#include <vector>

extern "C" {
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde_avx2.h> // AVX2 and lower
}

#include "opal.h"


// I define aliases for SSE intrinsics, so they can be used in code not depending on SSE generation.
// If available, AVX2 is used because it has two times bigger register, thus everything is two times faster.
// #ifdef __AVX2__
#if 1     //A.D.: now SIMDe will be taking care of the SIMD

const int SIMD_REG_SIZE = 256; //!< number of bits in register
typedef __m256i __mxxxi; //!< represents register containing integers
#define _mmxxx_load_si  _mm256_load_si256
#define _mmxxx_store_si _mm256_store_si256
#define _mmxxx_and_si   _mm256_and_si256
#define _mmxxx_testz_si _mm256_testz_si256

#define _mmxxx_adds_epi8  _mm256_adds_epi8
#define _mmxxx_subs_epi8  _mm256_subs_epi8
#define _mmxxx_min_epu8   _mm256_min_epu8
#define _mmxxx_min_epi8   _mm256_min_epi8
#define _mmxxx_max_epu8   _mm256_max_epu8
#define _mmxxx_max_epi8   _mm256_max_epi8
#define _mmxxx_set1_epi8  _mm256_set1_epi8
#define _mmxxx_cmpgt_epi8 _mm256_cmpgt_epi8

#define _mmxxx_adds_epi16  _mm256_adds_epi16
#define _mmxxx_subs_epi16  _mm256_subs_epi16
#define _mmxxx_min_epi16   _mm256_min_epi16
#define _mmxxx_max_epi16   _mm256_max_epi16
#define _mmxxx_set1_epi16  _mm256_set1_epi16
#define _mmxxx_cmpgt_epi16 _mm256_cmpgt_epi16

#define _mmxxx_add_epi32   _mm256_add_epi32
#define _mmxxx_sub_epi32   _mm256_sub_epi32
#define _mmxxx_min_epi32   _mm256_min_epi32
#define _mmxxx_max_epi32   _mm256_max_epi32
#define _mmxxx_set1_epi32  _mm256_set1_epi32
#define _mmxxx_cmpgt_epi32 _mm256_cmpgt_epi32

#else // SSE4.1

const int SIMD_REG_SIZE = 128;
typedef __m128i __mxxxi;
#define _mmxxx_load_si  _mm_load_si128
#define _mmxxx_store_si _mm_store_si128
#define _mmxxx_and_si   _mm_and_si128
#define _mmxxx_testz_si _mm_testz_si128

#define _mmxxx_adds_epi8  _mm_adds_epi8
#define _mmxxx_subs_epi8  _mm_subs_epi8
#define _mmxxx_min_epu8   _mm_min_epu8
#define _mmxxx_min_epi8   _mm_min_epi8
#define _mmxxx_max_epu8   _mm_max_epu8
#define _mmxxx_max_epi8   _mm_max_epi8
#define _mmxxx_set1_epi8  _mm_set1_epi8
#define _mmxxx_cmpgt_epi8 _mm_cmpgt_epi8

#define _mmxxx_adds_epi16  _mm_adds_epi16
#define _mmxxx_subs_epi16  _mm_subs_epi16
#define _mmxxx_min_epi16   _mm_min_epi16
#define _mmxxx_max_epi16   _mm_max_epi16
#define _mmxxx_set1_epi16  _mm_set1_epi16
#define _mmxxx_cmpgt_epi16 _mm_cmpgt_epi16

#define _mmxxx_add_epi32   _mm_add_epi32
#define _mmxxx_sub_epi32   _mm_sub_epi32
#define _mmxxx_min_epi32   _mm_min_epi32
#define _mmxxx_max_epi32   _mm_max_epi32
#define _mmxxx_set1_epi32  _mm_set1_epi32
#define _mmxxx_cmpgt_epi32 _mm_cmpgt_epi32

#endif


//------------------------------------ SIMD PARAMETERS ---------------------------------//
static inline int simdIsAllZeroes(const __mxxxi& a) {
    return _mmxxx_testz_si(a, a);
}

/**
 * Contains parameters and SIMD instructions specific for certain score type.
 */

template<typename T> struct SimdSW {};

template<>
struct SimdSW<char> {
    typedef char type; //!< Type that will be used for score
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(char)); //!< Number of sequences that can be done in parallel.
    static const bool satArthm = true; //!< True if saturation arithmetic is used, false otherwise.
    static const bool negRange = true; //!< True if it uses negative range for score representation, goes with saturation
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_adds_epi8(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_subs_epi8(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epu8(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epu8(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi8(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi8(a); }
};

template<>
struct SimdSW<short> {
    typedef short type;
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(short));
    static const bool satArthm = true;
    static const bool negRange = false;
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_adds_epi16(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_subs_epi16(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epi16(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epi16(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi16(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi16(a); }
};

template<>
struct SimdSW<int> {
    typedef int type;
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(int));
    static const bool satArthm = false;
    static const bool negRange = false;
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_add_epi32(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_sub_epi32(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epi32(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epi32(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi32(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi32(a); }
};
//--------------------------------------------------------------------------------------//


static bool loadNextSequence(int &nextDbSeqIdx, int dbLength, int &currDbSeqIdx, unsigned char * &currDbSeqPos,
                             int &currDbSeqLength, unsigned char ** db, int dbSeqLengths[], bool calculated[],
                             int &numEndedDbSeqs);


// For debugging
template<class SIMD>
void print_mmxxxi(__mxxxi mm) {
    typename SIMD::type unpacked[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
    _mmxxx_store_si((__mxxxi*)unpacked, mm);
    for (int i = 0; i < SIMD::numSeqs; i++)
        printf("%d ", unpacked[i]);
}

// This structure represents cell in dynamic programming matrix, but only with E and H values.
struct CellEH {
    __mxxxi H;
    __mxxxi E;
};

/**
 * @param stopOnOverflow  If true, function will stop when first overflow happens.
 *            If false, function will not stop but continue with next sequence.
 *
 */
template<class SIMD>
static int searchDatabaseSW_(unsigned char query[], int queryLength,
                             unsigned char** db, int dbLength, int dbSeqLengths[],
                             int gapOpen, int gapExt, int* scoreMatrix, int alphabetLength,
                             OpalSearchResult* results[], const int searchType, bool calculated[],
                             int overflowMethod) {

    const typename SIMD::type LOWER_BOUND = std::numeric_limits<typename SIMD::type>::min();
    const typename SIMD::type UPPER_BOUND = std::numeric_limits<typename SIMD::type>::max();

    bool overflowOccured = false;  // True if overflow was detected at least once.

    // ----------------------- CHECK ARGUMENTS -------------------------- //
    // Check if Q, R or scoreMatrix have values too big for used score type
    if (gapOpen < LOWER_BOUND || UPPER_BOUND < gapOpen || gapExt < LOWER_BOUND || UPPER_BOUND < gapExt) {
        return 1;
    }
    if (!SIMD::satArthm) {
        // These extra limits are enforced so overflow could be detected more efficiently
        if (gapOpen <= LOWER_BOUND/2 || UPPER_BOUND/2 <= gapOpen || gapExt <= LOWER_BOUND/2 || UPPER_BOUND/2 <= gapExt) {
            return OPAL_ERR_OVERFLOW;
        }
    }

    for (int r = 0; r < alphabetLength; r++)
        for (int c = 0; c < alphabetLength; c++) {
            int score = scoreMatrix[r * alphabetLength + c];
            if (score < LOWER_BOUND || UPPER_BOUND < score) {
                return OPAL_ERR_OVERFLOW;
            }
            if (!SIMD::satArthm) {
                if (score <= LOWER_BOUND/2 || UPPER_BOUND/2 <= score)
                    return OPAL_ERR_OVERFLOW;
            }
        }
    // ------------------------------------------------------------------ //


    // ------------------------ INITIALIZATION -------------------------- //
    __mxxxi zeroes = SIMD::set1(0);
    __mxxxi scoreZeroes; // 0 normally, but lower bound if using negative range
    if (SIMD::negRange) {
        scoreZeroes = SIMD::set1(LOWER_BOUND);
    } else {
        scoreZeroes = zeroes;
    }

    int numEndedDbSeqs = 0; // Number of sequences that ended
    int nextDbSeqIdx = 0;  // Index in db of next sequence that we will try to load.
    // Index in db. -1 for null sequence (null sequence == no sequence).
    int currDbSeqsIdxs[SIMD::numSeqs];
    // Pointer to current element for each current database sequence. 0 for null sequence.
    unsigned char* currDbSeqsPos[SIMD::numSeqs];
    // Lengths of remaining parts of sequences. -1 for null sequence.
    int currDbSeqsLengths[SIMD::numSeqs];

    // Needed in order to find the most early result, which is a nice condition to have.
    int currDbSeqsBestScore[SIMD::numSeqs];
    // Row index of best score for each current database sequence.
    int currDbSeqsBestScoreRow[SIMD::numSeqs];
    // Column index of best score for each current database sequence.
    int currDbSeqsBestScoreColumn[SIMD::numSeqs];

    // Profile query -> here we store preprocessed score data needed in core loop.
    // It is recalculated for each column.
    __mxxxi P[alphabetLength];

    // Load initial sequences
    for (int i = 0; i < SIMD::numSeqs; i++) {
        currDbSeqsBestScore[i] = LOWER_BOUND;
        loadNextSequence(nextDbSeqIdx, dbLength, currDbSeqsIdxs[i], currDbSeqsPos[i],
                         currDbSeqsLengths[i], db, dbSeqLengths, calculated, numEndedDbSeqs);
    }

    // Q is gap open penalty, R is gap ext penalty.
    __mxxxi Q = SIMD::set1(gapOpen);
    __mxxxi R = SIMD::set1(gapExt);

    int rowsWithImprovement[queryLength];  // Indexes of rows where one of sequences improved score.

    // Previous Hs, previous Es, previous F, all signed short.
    CellEH prevColumn[queryLength];  // Stores results of previous column in matrix.
    // Initialize all values to 0
    for (int i = 0; i < queryLength; i++) {
        prevColumn[i].H = prevColumn[i].E = scoreZeroes;
    }

    __mxxxi maxH = scoreZeroes;  // Best score in sequence
    // ------------------------------------------------------------------ //


    // For each column
    while (numEndedDbSeqs < dbLength) {
        // -------------------- CALCULATE QUERY PROFILE ------------------------- //
        // TODO: Rognes uses pshufb here, I don't know how/why?
        typename SIMD::type profileRow[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
        for (unsigned char letter = 0; letter < alphabetLength; letter++) {
            int* scoreMatrixRow = scoreMatrix + letter*alphabetLength;
            for (int i = 0; i < SIMD::numSeqs; i++) {
                unsigned char* dbSeqPos = currDbSeqsPos[i];
                if (dbSeqPos != 0)
                    profileRow[i] = (typename SIMD::type)scoreMatrixRow[*dbSeqPos];
            }
            P[letter] = _mmxxx_load_si((__mxxxi const*)profileRow);
        }
        // ---------------------------------------------------------------------- //

        // Previous cells: u - up, l - left, ul - up left
        __mxxxi uF, uH, ulH;
        uF = uH = ulH = scoreZeroes; // F[-1, c] = H[-1, c] = H[-1, c-1] = 0

        __mxxxi ofTest = scoreZeroes; // Used for detecting the overflow when not using saturated ar

        int rowsWithImprovementLength = 0;

        // ----------------------- CORE LOOP (ONE COLUMN) ----------------------- //
        for (int r = 0; r < queryLength; r++) { // For each cell in column
            // Calculate E = max(lH-Q, lE-R)
            __mxxxi E = SIMD::max(SIMD::sub(prevColumn[r].H, Q), SIMD::sub(prevColumn[r].E, R));

            // Calculate F = max(uH-Q, uF-R)
            __mxxxi F = SIMD::max(SIMD::sub(uH, Q), SIMD::sub(uF, R));

            // Calculate H
            __mxxxi H = SIMD::max(F, E);
            if (!SIMD::negRange) {
                // If not using negative range, then H could be negative at this moment so we need this
                H = SIMD::max(H, zeroes);
            }
            __mxxxi ulH_P = SIMD::add(ulH, P[query[r]]);
            // If using negative range: if ulH_P >= 0 then we have overflow

            H = SIMD::max(H, ulH_P);
            // If using negative range: H will always be negative, even if ulH_P overflowed

            // Save data needed for overflow detection. Not more then one condition will fire
            if (SIMD::negRange)
                ofTest = _mmxxx_and_si(ofTest, ulH_P);
            if (!SIMD::satArthm)
                ofTest = SIMD::min(ofTest, ulH_P);

            // If we need end location, remember row with best score.
            if (searchType != OPAL_SEARCH_SCORE) {
                // We remember rows that had max scores, in order to find out the row of best score.
                rowsWithImprovement[rowsWithImprovementLength] = r;
                // TODO(martin): simdIsAllZeroes seems to bring significant slowdown, but I could
                // not find way to avoid it.
                rowsWithImprovementLength += 1 - simdIsAllZeroes(SIMD::cmpgt(H, maxH));
            }

            maxH = SIMD::max(maxH, H); // update best score

            // Set uF, uH, ulH
            uF = F;
            uH = H;
            ulH = prevColumn[r].H;

            // Remeber values so they can be used in next column.
            prevColumn[r].E = E;
            prevColumn[r].H = H;

            // For saturated: score is biased everywhere, but just score: E, F, H
            // Also, all scores except ulH_P certainly have value < 0
        }
        // ---------------------------------------------------------------------- //

        for (int seqIdx = 0; seqIdx < SIMD::numSeqs; seqIdx++) {
            currDbSeqsLengths[seqIdx] -= currDbSeqsLengths[seqIdx] > 0;
        }

        typename SIMD::type unpackedMaxH[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
        _mmxxx_store_si((__mxxxi*)unpackedMaxH, maxH);

        // ------------------------ OVERFLOW DETECTION -------------------------- //
        bool overflowed[SIMD::numSeqs];
        if (!SIMD::satArthm) {
            // This check is based on following assumptions:
            //  - overflow wraps
            //  - Q, R and all scores from scoreMatrix are between LOWER_BOUND/2 and UPPER_BOUND/2 exclusive
            typename SIMD::type unpackedOfTest[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
            _mmxxx_store_si((__mxxxi*)unpackedOfTest, ofTest);
            for (int i = 0; i < SIMD::numSeqs; i++) {
                overflowed[i] = currDbSeqsPos[i] != 0 && unpackedOfTest[i] <= LOWER_BOUND / 2;
                if (overflowMethod == OPAL_OVERFLOW_BUCKETS && overflowed[i]) {
                    // In buckets method, we stop calculation when overflow is detected.
                    return OPAL_ERR_OVERFLOW;
                }
            }
        } else {
            if (SIMD::negRange) {
                // Since I use saturation, I check if minUlH_P was non negative
                typename SIMD::type unpackedOfTest[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
                _mmxxx_store_si((__mxxxi*)unpackedOfTest, ofTest);
                for (int i = 0; i < SIMD::numSeqs; i++) {
                    overflowed[i] = currDbSeqsPos[i] != 0 && unpackedOfTest[i] >= 0;
                    if (overflowMethod == OPAL_OVERFLOW_BUCKETS && overflowed[i]) {
                        return OPAL_ERR_OVERFLOW;
                    }
                }
            } else {
                // I check if upper bound is reached
                for (int i = 0; i < SIMD::numSeqs; i++) {
                    overflowed[i] = currDbSeqsPos[i] != 0 && unpackedMaxH[i] == UPPER_BOUND;
                    if (overflowMethod == OPAL_OVERFLOW_BUCKETS && overflowed[i]) {
                        return OPAL_ERR_OVERFLOW;
                    }
                }
            }
        }
        // ---------------------------------------------------------------------- //

        bool overflowDetected = false;  // True if overflow was detected in last column.
        bool sequenceEnded = false;  // True if a sequence ended in last column.
        for (int seqIdx = 0; seqIdx < SIMD::numSeqs; seqIdx++) {
            overflowDetected = overflowDetected || overflowed[seqIdx];
            sequenceEnded = sequenceEnded || (currDbSeqsLengths[seqIdx] == 0);
        }
        overflowOccured = overflowOccured || overflowDetected;


        // --------- Update end location of alignment ----------- //
        if (searchType != OPAL_SEARCH_SCORE) {
            for (int i = 0; i < rowsWithImprovementLength; i++) {
                int r = rowsWithImprovement[i];
                typename SIMD::type unpackedH[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
                _mmxxx_store_si((__mxxxi*)unpackedH, prevColumn[r].H);
                for (int j = 0; j < SIMD::numSeqs; j++) {
                    if (currDbSeqsPos[j] != 0 && !overflowed[j]) {  // If not null sequence or overflowedresult->endLocationQuery =
                        if (unpackedH[j] > currDbSeqsBestScore[j]) {
                            currDbSeqsBestScore[j] = unpackedH[j];
                            currDbSeqsBestScoreRow[j] = r;
                            currDbSeqsBestScoreColumn[j] = dbSeqLengths[currDbSeqsIdxs[j]]
                                - currDbSeqsLengths[j] - 1;
                        }
                    }
                }
            }
        }
        // ------------------------------------------------------ //


        // --------------------- CHECK AND HANDLE SEQUENCE END ------------------ //
        if (overflowDetected || sequenceEnded) { // If at least one sequence ended
            typename SIMD::type resetMask[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));

            for (int i = 0; i < SIMD::numSeqs; i++) {
                if (currDbSeqsPos[i] != 0) { // If not null sequence
                    if (overflowed[i] || currDbSeqsLengths[i] == 0) { // If sequence ended
                        numEndedDbSeqs++;
                        if (!overflowed[i]) {
                            // Save score and mark as calculated
                            calculated[currDbSeqsIdxs[i]] = true;
                            opalSearchResultSetScore(results[currDbSeqsIdxs[i]], unpackedMaxH[i]);
                            if (SIMD::negRange) {
                                results[currDbSeqsIdxs[i]]->score -= LOWER_BOUND;
                            }
                            if (searchType != OPAL_SEARCH_SCORE) {
                                results[currDbSeqsIdxs[i]]->endLocationQuery = currDbSeqsBestScoreRow[i];
                                results[currDbSeqsIdxs[i]]->endLocationTarget = currDbSeqsBestScoreColumn[i];
                            } else {
                                results[currDbSeqsIdxs[i]]->endLocationQuery = -1;
                                results[currDbSeqsIdxs[i]]->endLocationTarget = -1;
                            }
                        }
                        currDbSeqsBestScore[i] = LOWER_BOUND;
                        // Load next sequence
                        loadNextSequence(nextDbSeqIdx, dbLength, currDbSeqsIdxs[i], currDbSeqsPos[i],
                                         currDbSeqsLengths[i], db, dbSeqLengths, calculated, numEndedDbSeqs);
                        // If negative range, sets to LOWER_BOUND when used with saturated add and value < 0,
                        // otherwise sets to zero when used with and.
                        resetMask[i] = SIMD::negRange ? LOWER_BOUND : 0;
                    } else {
                        // Does not change anything when used with saturated add / and.
                        resetMask[i] = SIMD::negRange ? 0 : -1;
                        if (currDbSeqsPos[i] != 0)
                            currDbSeqsPos[i]++; // If not new and not null, move for one element
                    }
                }
            }
            // Reset previous columns (Es and Hs) and maxH
            __mxxxi resetMaskPacked = _mmxxx_load_si((__mxxxi const*)resetMask);
            if (SIMD::negRange) {
                for (int i = 0; i < queryLength; i++) {
                    prevColumn[i].H = SIMD::add(prevColumn[i].H, resetMaskPacked);
                    prevColumn[i].E = SIMD::add(prevColumn[i].E, resetMaskPacked);
                }
                maxH = SIMD::add(maxH, resetMaskPacked);
            } else {
                for (int i = 0; i < queryLength; i++) {
                    prevColumn[i].H = _mmxxx_and_si(prevColumn[i].H, resetMaskPacked);
                    prevColumn[i].E = _mmxxx_and_si(prevColumn[i].E, resetMaskPacked);
                }
                maxH = _mmxxx_and_si(maxH, resetMaskPacked);
            }
        } else { // If no sequences ended
            // Move for one element in all sequences
            for (int i = 0; i < SIMD::numSeqs; i++)
                currDbSeqsPos[i] += (currDbSeqsPos[i] != 0);
        }
        // ---------------------------------------------------------------------- //
    }

    if (overflowOccured) {
        return OPAL_ERR_OVERFLOW;
    }
    return 0;
}

static inline bool loadNextSequence(int &nextDbSeqIdx, int dbLength, int &currDbSeqIdx, unsigned char* &currDbSeqPos,
                                    int &currDbSeqLength, unsigned char** db, int dbSeqLengths[], bool calculated[],
                                    int &numEndedDbSeqs) {
    while (nextDbSeqIdx < dbLength && calculated[nextDbSeqIdx]) {
        nextDbSeqIdx++;
        numEndedDbSeqs++;
    }
    if (nextDbSeqIdx < dbLength) { // If there is sequence to load
        currDbSeqIdx = nextDbSeqIdx;
        currDbSeqPos = db[nextDbSeqIdx];
        currDbSeqLength = dbSeqLengths[nextDbSeqIdx];
        nextDbSeqIdx++;
        return true;
    } else { // If there are no more sequences to load, load "null" sequence
        currDbSeqIdx = currDbSeqLength = -1; // Set to -1 if there are no more sequences
        currDbSeqPos = 0;
        return false;
    }
}

/**
 * @param If skip[i] is true, result for sequence #i will not be calculated.
 *     If skip is NULL, all results will be calculated.
 */
static int searchDatabaseSW(unsigned char query[], int queryLength,
                            unsigned char** db, int dbLength, int dbSeqLengths[],
                            int gapOpen, int gapExt, int* scoreMatrix, int alphabetLength,
                            OpalSearchResult* results[], const int searchType, bool skip[], int overflowMethod) {
    int resultCode = 0;
    // Do buckets only if using buckets overflow method.
    const int chunkSize = overflowMethod == OPAL_OVERFLOW_BUCKETS ? 1024 : dbLength;
    bool* calculated = new bool[chunkSize];
    for (int startIdx = 0; startIdx < dbLength; startIdx += chunkSize) {
        unsigned char** db_ = db + startIdx;
        int* dbSeqLengths_ = dbSeqLengths + startIdx;
        OpalSearchResult** results_ = results + startIdx;
        int dbLength_ = startIdx + chunkSize >= dbLength ? dbLength - startIdx : chunkSize;
        for (int i = 0; i < dbLength_; i++) {
            calculated[i] = skip ? skip[i] : false;
        }
        resultCode = searchDatabaseSW_< SimdSW<char> >(
            query, queryLength, db_, dbLength_, dbSeqLengths_,
            gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
            searchType, calculated, overflowMethod);
        if (resultCode == OPAL_ERR_OVERFLOW) {
            resultCode = searchDatabaseSW_< SimdSW<short> >(
                query, queryLength, db_, dbLength_, dbSeqLengths_,
                gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
                searchType, calculated, overflowMethod);
            if (resultCode == OPAL_ERR_OVERFLOW) {
                resultCode = searchDatabaseSW_< SimdSW<int> >(
                    query, queryLength,
                    db_, dbLength_, dbSeqLengths_,
                    gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
                    searchType, calculated, overflowMethod);
                if (resultCode != 0)
                    break;
            }
        }
    }

    delete[] calculated;
    return resultCode;
}









//------------------------------------ SIMD PARAMETERS ---------------------------------//
/**
 * Contains parameters and SIMD instructions specific for certain score type.
 */
template<typename T> class Simd {};

template<>
struct Simd<char> {
    typedef char type; //!< Type that will be used for score
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(char)); //!< Number of sequences that can be done in parallel.
    static const bool satArthm = true; //!< True if saturation arithmetic is used, false otherwise.
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_adds_epi8(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_subs_epi8(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epi8(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epi8(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi8(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi8(a); }
};

template<>
struct Simd<short> {
    typedef short type;
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(short));
    static const bool satArthm = true;
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_adds_epi16(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_subs_epi16(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epi16(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epi16(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi16(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi16(a); }
};

template<>
struct Simd<int> {
    typedef int type;
    static const int numSeqs = SIMD_REG_SIZE / (8 * sizeof(int));
    static const bool satArthm = false;
    static inline __mxxxi add(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_add_epi32(a, b); }
    static inline __mxxxi sub(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_sub_epi32(a, b); }
    static inline __mxxxi min(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_min_epi32(a, b); }
    static inline __mxxxi max(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_max_epi32(a, b); }
    static inline __mxxxi cmpgt(const __mxxxi& a, const __mxxxi& b) { return _mmxxx_cmpgt_epi32(a, b); }
    static inline __mxxxi set1(int a) { return _mmxxx_set1_epi32(a); }
};
//--------------------------------------------------------------------------------------//




template<class SIMD, int MODE>
static int searchDatabase_(unsigned char query[], int queryLength,
                           unsigned char** db, int dbLength, int dbSeqLengths[],
                           int gapOpen, int gapExt, int* scoreMatrix, int alphabetLength,
                           OpalSearchResult* results[], const int searchType, bool calculated[], int overflowMethod) {

    static const typename SIMD::type LOWER_BOUND = std::numeric_limits<typename SIMD::type>::min();
    static const typename SIMD::type UPPER_BOUND = std::numeric_limits<typename SIMD::type>::max();
    // Used to represent -inf. Must be larger then lower bound to avoid overflow.
    // TOOD(martin): Is it enough to use gapExt here? Should I use gapOpen instead?
    static const typename SIMD::type LOWER_SCORE_BOUND = LOWER_BOUND + gapExt;

    bool overflowOccured = false;  // True if oveflow was detected at least once.

    // ----------------------- CHECK ARGUMENTS -------------------------- //
    // Check if Q, R or scoreMatrix have values too big for used score type
    if (gapOpen < LOWER_BOUND || UPPER_BOUND < gapOpen || gapExt < LOWER_BOUND || UPPER_BOUND < gapExt) {
        return 1;
    }
    if (!SIMD::satArthm) {
        // These extra limits are enforced so overflow could be detected more efficiently
        if (gapOpen <= LOWER_BOUND/2 || UPPER_BOUND/2 <= gapOpen || gapExt <= LOWER_BOUND/2 || UPPER_BOUND/2 <= gapExt) {
            return 1;
        }
    }

    for (int r = 0; r < alphabetLength; r++)
        for (int c = 0; c < alphabetLength; c++) {
            int score = scoreMatrix[r * alphabetLength + c];
            if (score < LOWER_BOUND || UPPER_BOUND < score) {
                return 1;
            }
            if (!SIMD::satArthm) {
                if (score <= LOWER_BOUND/2 || UPPER_BOUND/2 <= score)
                    return 1;
            }
        }

    // TODO: If not saturated arithmetic check if inital values of H (-gapOpen - i * gapExt) cause overflow
    // ------------------------------------------------------------------ //


    // ------------------------ INITIALIZATION -------------------------- //
    const __mxxxi ZERO_SIMD = SIMD::set1(0);
    const __mxxxi LOWER_BOUND_SIMD = SIMD::set1(LOWER_BOUND);
    const __mxxxi LOWER_SCORE_BOUND_SIMD = SIMD::set1(LOWER_SCORE_BOUND);

    int nextDbSeqIdx = 0; // index in db
    int currDbSeqsIdxs[SIMD::numSeqs]; // index in db
    unsigned char* currDbSeqsPos[SIMD::numSeqs]; // current element for each current database sequence
    int currDbSeqsLengths[SIMD::numSeqs];
    bool justLoaded[SIMD::numSeqs] = {0}; // True if sequence was just loaded into channel
    bool seqJustLoaded = false; // True if at least one sequence was just loaded into channel
    int shortestDbSeqLength = -1;  // length of shortest sequence among current database sequences
    int numEndedDbSeqs = 0; // Number of sequences that ended

    // Row index of best score for each current database sequence. Used for HW and OV.
    int currDbSeqsBestScoreColumn[SIMD::numSeqs];

    // Load initial sequences
    for (int i = 0; i < SIMD::numSeqs; i++)
        if (loadNextSequence(nextDbSeqIdx, dbLength, currDbSeqsIdxs[i], currDbSeqsPos[i],
                             currDbSeqsLengths[i], db, dbSeqLengths, calculated, numEndedDbSeqs)) {
            justLoaded[i] = seqJustLoaded = true;
            // Update shortest sequence length if new sequence was loaded
            if (shortestDbSeqLength == -1 || currDbSeqsLengths[i] < shortestDbSeqLength)
                shortestDbSeqLength = currDbSeqsLengths[i];
        }

    // Q is gap open penalty, R is gap ext penalty.
    const __mxxxi Q = SIMD::set1(gapOpen);
    const __mxxxi R = SIMD::set1(gapExt);

    // Previous H column (array), previous E column (array), previous F, all signed short
    __mxxxi prevHs[queryLength];
    __mxxxi prevEs[queryLength];
    // Initialize all values
    for (int r = 0; r < queryLength; r++) {
        if (MODE == OPAL_MODE_OV)
            prevHs[r] = ZERO_SIMD;
        else { // - Q - r * R
            if (r == 0)
                prevHs[0] = SIMD::sub(ZERO_SIMD, Q);
            else
                prevHs[r] = SIMD::sub(prevHs[r-1], R);
        }

        prevEs[r] = LOWER_SCORE_BOUND_SIMD;
    }

    // u - up, ul - up left
    __mxxxi uH, ulH;
    if (MODE == OPAL_MODE_NW) {
        ulH = ZERO_SIMD;
        uH = SIMD::sub(R, Q); // -Q + R
    }

    __mxxxi maxLastRowH = LOWER_BOUND_SIMD; // Keeps track of maximum H in last row
    // ------------------------------------------------------------------ //


    int columnsSinceLastSeqEnd = 0;
    // For each column
    while (numEndedDbSeqs < dbLength) {
        // -------------------- CALCULATE QUERY PROFILE ------------------------- //
        // TODO: Rognes uses pshufb here, I don't know how/why?
        __mxxxi P[alphabetLength];
        typename SIMD::type profileRow[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
        for (unsigned char letter = 0; letter < alphabetLength; letter++) {
            int* scoreMatrixRow = scoreMatrix + letter*alphabetLength;
            for (int i = 0; i < SIMD::numSeqs; i++) {
                unsigned char* dbSeqPos = currDbSeqsPos[i];
                if (dbSeqPos != 0)
                    profileRow[i] = (typename SIMD::type)scoreMatrixRow[*dbSeqPos];
            }
            P[letter] = _mmxxx_load_si((__mxxxi const*)profileRow);
        }
        // ---------------------------------------------------------------------- //

        // u - up
        __mxxxi uF = LOWER_SCORE_BOUND_SIMD;

        // Database sequence has fixed start and end only in NW
        if (MODE == OPAL_MODE_NW) {
            if (seqJustLoaded) {
                typename SIMD::type resetMask[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
                for (int i = 0; i < SIMD::numSeqs; i++)
                    resetMask[i] = justLoaded[i] ?  0 : -1;
                const __mxxxi resetMaskPacked = _mmxxx_load_si((__mxxxi const*)resetMask);
                ulH = _mmxxx_and_si(uH, resetMaskPacked);
            } else {
                ulH = uH;
            }

            uH = SIMD::sub(uH, R); // uH is -Q - c*R
            // NOTE: Setup of ulH and uH for first column is done when sequence is loaded.
        } else {
            uH = ulH = ZERO_SIMD;
        }

        __mxxxi minE, minF;
        minE = minF = SIMD::set1(UPPER_BOUND);
        __mxxxi maxH = LOWER_BOUND_SIMD; // Max H in this column
        __mxxxi H;

        __mxxxi firstRow_uH, firstRow_ulH; // Values of uH and ulH from first row of column

        if (MODE == OPAL_MODE_NW) {
            firstRow_uH = uH;
            firstRow_ulH = ulH;
        }

        __mxxxi prevMaxLastRowH = maxLastRowH;
        // ----------------------- CORE LOOP (ONE COLUMN) ----------------------- //
        for (int r = 0; r < queryLength; r++) { // For each cell in column
            // Calculate E = max(lH-Q, lE-R)
            __mxxxi E = SIMD::max(SIMD::sub(prevHs[r], Q), SIMD::sub(prevEs[r], R)); // E could overflow

            // Calculate F = max(uH-Q, uF-R)
            __mxxxi F = SIMD::max(SIMD::sub(uH, Q), SIMD::sub(uF, R)); // F could overflow
            minF = SIMD::min(minF, F); // For overflow detection

            // Calculate H
            H = SIMD::max(F, E);
            __mxxxi ulH_P = SIMD::add(ulH, P[query[r]]);
            H = SIMD::max(H, ulH_P); // H could overflow

            maxH = SIMD::max(maxH, H); // update best score in column

            // Set uF, uH, ulH
            uF = F;
            uH = H;
            ulH = prevHs[r];

            // Update prevHs, prevEs in advance for next column
            prevEs[r] = E;
            prevHs[r] = H;
        }
        // ---------------------------------------------------------------------- //

        maxLastRowH = SIMD::max(maxLastRowH, H);

        if (MODE == OPAL_MODE_NW) {
            uH = firstRow_uH;
            ulH = firstRow_ulH;
        }

        // Find minE, which should be checked with minE == LOWER_BOUND for overflow
        for (int r = 0; r < queryLength; r++)
            minE = SIMD::min(minE, prevEs[r]);

        columnsSinceLastSeqEnd++;

        typename SIMD::type unpackedMaxH[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
        _mmxxx_store_si((__mxxxi*)unpackedMaxH, maxH);

        // ------------------------ OVERFLOW DETECTION -------------------------- //
        bool overflowDetected = false;  // True if overflow was detected for this column.
        bool overflowed[SIMD::numSeqs];
        if (!SIMD::satArthm) {
            /*           // This check is based on following assumptions:
            //  - overflow wraps
            //  - Q, R and all scores from scoreMatrix are between LOWER_BOUND/2 and UPPER_BOUND/2 exclusive
            typename SIMD::type* unpackedOfTest = (typename SIMD::type *)&ofTest;
            for (int i = 0; i < SIMD::numSeqs; i++)
                if (currDbSeqsPos[i] != 0 && unpackedOfTest[i] <= LOWER_BOUND/2)
                return 1;*/
        } else {
            // There is overflow if minE == LOWER_BOUND or minF == LOWER_BOUND or maxH == UPPER_BOUND
            __mxxxi minEF = SIMD::min(minE, minF);
            typename SIMD::type unpackedMinEF[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
            _mmxxx_store_si((__mxxxi*)unpackedMinEF, minEF);
            for (int i = 0; i < SIMD::numSeqs; i++) {
                overflowed[i] = currDbSeqsPos[i] != 0 && (unpackedMinEF[i] == LOWER_BOUND
                                                          || unpackedMaxH[i] == UPPER_BOUND);
                if (overflowMethod == OPAL_OVERFLOW_BUCKETS && overflowed[i]) {
                    // In buckets method, we stop calculation when overflow is detected.
                    return OPAL_ERR_OVERFLOW;
                }
            }
        }
        for (int i = 0; i < SIMD::numSeqs; i++) {
            overflowDetected = overflowDetected || overflowed[i];
        }
        overflowOccured = overflowOccured || overflowDetected;
        // ---------------------------------------------------------------------- //

        // ------------------ Store end location of best score ------------------ //
        if (searchType != OPAL_SEARCH_SCORE && (MODE == OPAL_MODE_HW || MODE == OPAL_MODE_OV)) {
            // Determine the column of best score.
            __mxxxi greater = SIMD::cmpgt(maxLastRowH, prevMaxLastRowH);
            typename SIMD::type unpackedGreater[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
            _mmxxx_store_si((__mxxxi*)unpackedGreater, greater);
            for (int i = 0; i < SIMD::numSeqs; i++) {
                if (currDbSeqsPos[i] != 0 && !overflowed[i]) {  // If not null sequence or overflowed
                    if (unpackedGreater[i] != 0) {
                        currDbSeqsBestScoreColumn[i] = dbSeqLengths[currDbSeqsIdxs[i]] - currDbSeqsLengths[i]
                            + columnsSinceLastSeqEnd - 1;
                    }
                }
            }
        }
        // ---------------------------------------------------------------------- //

        seqJustLoaded = false;
        // --------------------- CHECK AND HANDLE SEQUENCE END ------------------ //
        if (overflowDetected || shortestDbSeqLength == columnsSinceLastSeqEnd) { // If at least one sequence ended
            shortestDbSeqLength = -1;

            // Calculate best scores
            __mxxxi bestScore;
            if (MODE == OPAL_MODE_OV)
                bestScore = SIMD::max(maxH, maxLastRowH); // Maximum of last row and column
            if (MODE == OPAL_MODE_HW)
                bestScore = maxLastRowH;
            if (MODE == OPAL_MODE_NW)
                bestScore = H;
            typename SIMD::type unpackedBestScore[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
            _mmxxx_store_si((__mxxxi*)unpackedBestScore, bestScore);

            for (int i = 0; i < SIMD::numSeqs; i++) {
                if (currDbSeqsPos[i] != 0) { // If not null sequence
                    justLoaded[i] = false;
                    currDbSeqsLengths[i] -= columnsSinceLastSeqEnd;

                    if (overflowed[i] || currDbSeqsLengths[i] == 0) { // If sequence ended
                        numEndedDbSeqs++;
                        if (!overflowed[i]) {
                            // Store result.
                            int dbSeqIdx = currDbSeqsIdxs[i];
                            OpalSearchResult *result = results[dbSeqIdx];
                            calculated[dbSeqIdx] = true;
                            // Set score.
                            opalSearchResultSetScore(result, unpackedBestScore[i]);
                            // Set end location.
                            if (searchType == OPAL_SEARCH_SCORE) {
                                result->endLocationQuery = -1;
                                result->endLocationTarget = -1;
                            } else {
                                if (MODE == OPAL_MODE_NW) {
                                    result->endLocationQuery = queryLength - 1;
                                    result->endLocationTarget = dbSeqLengths[dbSeqIdx] - 1;
                                }
                                if (MODE == OPAL_MODE_HW) {
                                    result->endLocationQuery = queryLength - 1;
                                    result->endLocationTarget = currDbSeqsBestScoreColumn[i];
                                }
                                if (MODE == OPAL_MODE_OV) {
                                    // This unpacking will repeat unnecessarily if there are multiple sequences
                                    // ending at the same time, however that will happen in very rare occasions.
                                    // TODO(martin): always unpack only once.
                                    typename SIMD::type unpackedPrevMaxLastRowH[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
                                    _mmxxx_store_si((__mxxxi*)unpackedPrevMaxLastRowH, prevMaxLastRowH);
                                    typename SIMD::type maxScore = unpackedPrevMaxLastRowH[i];

                                    // If last column contains best score, determine row.
                                    if (unpackedMaxH[i] > maxScore) {
                                        result->endLocationTarget = dbSeqLengths[dbSeqIdx] - 1;
                                        for (int r = 0; r < queryLength; r++) {
                                            typename SIMD::type unpackedPrevH[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
                                            _mmxxx_store_si((__mxxxi*)unpackedPrevH, prevHs[r]);
                                            if (unpackedPrevH[i] > maxScore) {
                                                result->endLocationQuery = r;
                                                maxScore = unpackedPrevH[i];
                                            }
                                        }
                                    } else {
                                        result->endLocationTarget = currDbSeqsBestScoreColumn[i];
                                        result->endLocationQuery = queryLength - 1;
                                    }
                                }
                            }
                        }
                        // Load next sequence
                        if (loadNextSequence(nextDbSeqIdx, dbLength, currDbSeqsIdxs[i], currDbSeqsPos[i],
                                             currDbSeqsLengths[i], db, dbSeqLengths, calculated, numEndedDbSeqs)) {
                            justLoaded[i] = seqJustLoaded = true;
                        }
                    } else {
                        if (currDbSeqsPos[i] != 0)
                            currDbSeqsPos[i]++; // If not new and not null, move for one element
                    }
                    // Update shortest sequence length if sequence is not null
                    if (currDbSeqsPos[i] != 0 && (shortestDbSeqLength == -1 || currDbSeqsLengths[i] < shortestDbSeqLength))
                        shortestDbSeqLength = currDbSeqsLengths[i];
                }
            }
            //------------ Reset prevEs, prevHs, maxLastRowH(, ulH and uH) ------------//
            typename SIMD::type resetMask[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8)));
            typename SIMD::type setMask[SIMD::numSeqs] __attribute__((aligned(SIMD_REG_SIZE / 8))); // inverse of resetMask
            for (int i = 0; i < SIMD::numSeqs; i++) {
                resetMask[i] = justLoaded[i] ?  0 : -1;
                setMask[i]   = justLoaded[i] ? -1 :  0;
            }
            const __mxxxi resetMaskPacked = _mmxxx_load_si((__mxxxi const*)resetMask);
            const __mxxxi setMaskPacked = _mmxxx_load_si((__mxxxi const*)setMask);

            // Set prevEs ended channels to LOWER_SCORE_BOUND
            const __mxxxi maskedLowerScoreBoundSimd = _mmxxx_and_si(setMaskPacked, LOWER_SCORE_BOUND_SIMD);
            for (int r = 0; r < queryLength; r++) {
                prevEs[r] = _mmxxx_and_si(prevEs[r], resetMaskPacked);
                prevEs[r] = SIMD::add(prevEs[r], maskedLowerScoreBoundSimd);
            }

            // Set prevHs
            for (int r = 0; r < queryLength; r++) {
                prevHs[r] = _mmxxx_and_si(prevHs[r], resetMaskPacked);
                if (MODE != OPAL_MODE_OV) {
                    if (r == 0) {
                        prevHs[0] = SIMD::sub(prevHs[0], _mmxxx_and_si(setMaskPacked, Q));
                    } else {
                        prevHs[r] = SIMD::add(prevHs[r], _mmxxx_and_si(setMaskPacked, SIMD::sub(prevHs[r-1], R)));
                    }
                }
            }

            // Set ulH and uH if NW
            if (MODE == OPAL_MODE_NW) {
                ulH = _mmxxx_and_si(ulH, resetMaskPacked); // to 0
                // Set uH channels to -Q + R
                uH = _mmxxx_and_si(uH, resetMaskPacked);
                uH = SIMD::add(uH, _mmxxx_and_si(setMaskPacked, SIMD::sub(R, Q)));
            }

            // Set maxLastRow ended channels to LOWER_BOUND
            maxLastRowH = _mmxxx_and_si(maxLastRowH, resetMaskPacked);
            maxLastRowH = SIMD::add(maxLastRowH, _mmxxx_and_si(setMaskPacked, LOWER_BOUND_SIMD));
            //-------------------------------------------------------//

            columnsSinceLastSeqEnd = 0;
        } else { // If no sequences ended
            // Move for one element in all sequences
            for (int i = 0; i < SIMD::numSeqs; i++)
                if (currDbSeqsPos[i] != 0)
                    currDbSeqsPos[i]++;
        }
        // ---------------------------------------------------------------------- //
    }

    if (overflowOccured) {
        return OPAL_ERR_OVERFLOW;
    }
    return 0;
}

/**
 * @param [in] skip  If skip[i] is true, result for sequence #i will not be calculated.
 *     If skip is NULL, all results will be calculated.
 */
template <int MODE>
static int searchDatabase(unsigned char query[], int queryLength,
                          unsigned char** db, int dbLength, int dbSeqLengths[],
                          int gapOpen, int gapExt, int* scoreMatrix, int alphabetLength,
                          OpalSearchResult* results[], const int searchType, bool skip[], int overflowMethod) {
    int resultCode = 0;
    // Do buckets only if using buckets overflow method.
    const int chunkSize = overflowMethod == OPAL_OVERFLOW_BUCKETS ? 1024 : dbLength;
    bool* calculated = new bool[chunkSize];
    for (int startIdx = 0; startIdx < dbLength; startIdx += chunkSize) {
        unsigned char** db_ = db + startIdx;
        int* dbSeqLengths_ = dbSeqLengths + startIdx;
        OpalSearchResult** results_ = results + startIdx;
        int dbLength_ = startIdx + chunkSize >= dbLength ? dbLength - startIdx : chunkSize;
        for (int i = 0; i < dbLength_; i++) {
            calculated[i] = skip ? skip[i] : false;
        }
        resultCode = searchDatabase_< Simd<char>, MODE >
            (query, queryLength, db_, dbLength_, dbSeqLengths_,
             gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
             searchType, calculated, overflowMethod);
        if (resultCode == OPAL_ERR_OVERFLOW) {
            resultCode = searchDatabase_< Simd<short>, MODE >
                (query, queryLength, db_, dbLength_, dbSeqLengths_,
                 gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
                 searchType, calculated, overflowMethod);
            if (resultCode == OPAL_ERR_OVERFLOW) {
                resultCode = searchDatabase_< Simd<int>, MODE >
                    (query, queryLength, db_, dbLength_, dbSeqLengths_,
                     gapOpen, gapExt, scoreMatrix, alphabetLength, results_,
                     searchType, calculated, overflowMethod);
                if (resultCode != 0)
                    break; // TODO: this does not make much sense because of buckets, improve it.
            }
        }
    }
    delete[] calculated;
    return resultCode;
}


/**
 * @return Zero-based index of max element in array. If there are multiple elements with
 *     same value, last is returned.
 */
template <typename T>
static inline T arrayMax(T array[], int length) {
    assert(length > 0);
    int maxElementIdx = 0;
    for (int i = 1; i < length; i++) {
        if (array[i] > array[maxElementIdx]) {
            maxElementIdx = i;
        }
    }
    return array[maxElementIdx];
}

/**
 * @param length  Length of gap, must be non-negative.
 * @param gapOpen  Penalty for gap opening (non-negative).
 * @param gapExt  Penalty for gap extension (non-negative).
 * @return Gap penalty, as non-negative number.
 */
static int gapPenalty(int length, int gapOpen, int gapExt) {
    if (length > 0) {
        return gapOpen + gapExt * (length - 1);
    } else {
        return 0;
    }
}

/**
 * Calculates bottom border of band for OV mode stop conditions.
 */
static int calculateBottomBandBorderOV(int k, int Q, int T, int Go, int Ge, int M) {
    int border = 0;

    // for d <= Q - T:
    border = std::max(border, std::min(Q - T, -1 * (k + Go - Ge - M * T) / Ge));

    // for d > Q - T:
    int borderCandidate = -1 * (k - M * Q + Go - Ge) / (Ge + M);
    if (borderCandidate > Q - T) {
        border = std::max(border, borderCandidate);
    }

    return std::min(border, Q - 1);
}

static int calculateTopBandBorderHW(int k, int Q, int T, int Go, int Ge, int M) {
    int border = 0;

    // for d <= T - Q;
    border = std::max(border, std::min(T - Q, -1 * (k - M * Q + Go) / Ge + 1));

    // for d > T - Q:
    int borderCandidate = -1 * (k - T * M + 2 * Go + Ge * (Q - T - 2)) / (2 * Ge + M);
    if (borderCandidate > T - Q) {
        border = std::max(border, borderCandidate);
    }

    return std::min(border, T - 1);
}

static int calculateBottomBandBorderHW(int k, int Q, int T, int Go, int Ge, int M) {
    int border = 0;

    // for d >= Q - T:
    int borderCandidate = -1 * (k + Go - Ge - Q * M) / (Ge + M);
    if (borderCandidate >= Q - T) {
        border = std::max(border, borderCandidate);
    }

    // for d < Q - T:
    if (-2 * Go - Ge * (Q - T - 2) + M * T >= k) {
        border = std::max(border, Q - T - 1);
    }

    return std::min(border, Q - 1);
}

static int calculateBottomBandBorderNW(int k, int Q, int T, int Go, int Ge, int M) {
    int border = 0;

    // for d > Q - T:
    int borderCandidate = -1 * (k + 2 * Go - M * Q + Ge * (T - Q - 2)) / (2 * Ge + M);
    if (borderCandidate > Q - T) {
        border = std::max(border, borderCandidate);
    }

    // for d = Q - T:
    if (Q - T <= -1 * (k + Go - M * T - Ge) / Ge) {
        border = std::max(border, Q - T);
    }

    // for d < Q - T:
    if (-2 * Go - Ge * (Q - T - 2) + M * T >= k) {
        border = std::max(border, Q - T - 1);
    }

    return std::min(border, Q - 1);
}

/**
 * Calculates top and bottom diagonal of band that contains all cells that could be part of any
 * solution that gives score not smaller than k. This means that if we are interested only in solutions
 * that give score that is not smaller than k, it is enough to calculate only cells inside that band.
 * Always starts from top left corner (like NW does) no matter which mode is specified,
 * and stops with regard to stop conditions of specified mode.
 * @param k  We are interested only in scores thar are not smaller than k.
 * @param mode  Alignment mode, it is used for stop conditions.
 * @param Q  queryLength
 * @param T  targetLength
 * @param Go  gapOpen -> non negative penalty for opening of gap.
 * @param Ge  gapExt -> non negative penalty for extension of gap.
 * @param M  Max score from score matrix.
 * @return Pair where first is index of bottom diagonal(border), and second is index of top diagonal(border).
 *     Therefore, first value will be in [0, Q - 1], while second will be in [0, T - 1].
 *     Band spans between top and bottom diagonal, including them.
 *     Main diagonal has index 0, diagonals above it and below both start with index 1.
 *     If there is no band, (-1, -1) is returned.
 *
 *     Example of matrix where Q is 3 and T is 6, each cell is marked with index of diagonal it lies on:
 *
 *          012345                                                         xxx---
 *          101234   , if bottom border is 1 and top is 2, we have band:   xxxx--
 *          210123                                                         -xxxx-
 */
static std::pair<int, int> calculateBandBorders(int k, int mode, int Q, int T, int Go, int Ge, int M) {
    if (mode == OPAL_MODE_OV || mode == OPAL_MODE_SW) {
        // Bands for OV and SW have same conditions, so they are calculated in same way.
        if (M * std::min(Q, T) >= k) {  // Determine if band exists at all.
            // Conditions for top and bottom band are symmetric, so logic for bottom is reused for top.
            return std::make_pair(calculateBottomBandBorderOV(k, Q, T, Go, Ge, M),
                                  calculateBottomBandBorderOV(k, T, Q, Go, Ge, M));
        } else {
            return std::make_pair(-1, -1);
        }
    } else if (mode == OPAL_MODE_HW) {
        if (M * std::min(Q, T) - gapPenalty(Q - std::min(Q, T), Go, Ge) >= k) {
            return std::make_pair(calculateBottomBandBorderHW(k, Q, T, Go, Ge, M),
                                  calculateTopBandBorderHW(k, Q, T, Go, Ge, M));
        } else {
            return std::make_pair(-1, -1);
        }
    } else if (mode == OPAL_MODE_NW) {
        if (M * std::min(Q, T) - gapPenalty(std::abs(Q - T), Go, Ge) >= k) {
            // Conditions for top and bottom band are symmetric, so logic for bottom is reused for top.
            return std::make_pair(calculateBottomBandBorderNW(k, Q, T, Go, Ge, M),
                                  calculateBottomBandBorderNW(k, T, Q, Go, Ge, M));
        } else {
            return std::make_pair(-1, -1);
        }
    } else {
        assert(false);  // Invalid alignment mode.
    }
}



/**
 * Returns new sequence that is reverse of given sequence.
 */
static inline unsigned char* createReverseCopy(const unsigned char* seq, int length) {
    unsigned char* rSeq = (unsigned char*) malloc(length * sizeof(unsigned char));
    for (int i = 0; i < length; i++) {
        rSeq[i] = seq[length - i - 1];
    }
    return rSeq;
}

template <class T>
static inline void revertArray(T array[], int length) {
    for (int i = 0; i < length / 2; i++) {
        T tmp = array[i];
        array[i] = array[length - 1 - i];
        array[length - 1 - i] = tmp;
    }
}

// Here I store scores for one cell in score matrix.
class Cell {
public:
    int H, E, F;
    enum class Field {
        H, E, F
    };
};

/**
 * Finds alignment of two sequences, if we know scoreLimit.
 * First alignment that has score greater or equal then scoreLimit will be returned.
 * If there is no such alignment, behavior will be unexpected.
 * Always starts from top left corner (like NW does) no matter which mode is specified,
 * and stops with regard to stop conditions of specified mode.
 * For example, for HW it will stop on last row, and for SW it will stop anywhere.
 * Returns score, start location (which is always (0, 0)), end location and alignment.
 * @param [in] query
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] gapOpen
 * @param [in] gapExt
 * @param [in] scoreMatrix
 * @param [in] alphabetLength
 * @param [in] scoreLimit  First alignment with score greater/equal than scoreLimit is returned.
 *     If there is no such score, behavior is undefined.
 *     TODO(martin): make this function also work when max score is smaller then scoreLimit.
 * @param [out] result  Pointer to already allocated object is expected here.
 *     Score, start location, end location and alignment will be set.
 *     Do not forget to free() alignment!
 * @param [in] mode  Mode whose stop conditions will be used when finding alignment.
 */
static void findAlignment(
    const unsigned char query[], const int queryLength, const unsigned char target[], const int targetLength,
    const int gapOpen, const int gapExt, const int* scoreMatrix, const int alphabetLength,
    const int scoreLimit, OpalSearchResult* result, const int mode) {
    /*
    printf("Query: ");
    for (int i = 0; i < queryLength; i++) {
        printf("%d ", query[i]);
    }
    printf("\n");
    printf("Target: ");
    for (int i = 0; i < targetLength; i++) {
        printf("%d ", target[i]);
    }
    printf("\n");
    */

    //printf("lengths: %d %d\n", queryLength, targetLength);

    std::pair<int, int> bandBorders = calculateBandBorders(
        scoreLimit, mode, queryLength, targetLength, gapOpen, gapExt,
        arrayMax(scoreMatrix, alphabetLength * alphabetLength));

    assert(bandBorders.first >= 0 && bandBorders.first < queryLength);
    assert(bandBorders.second >= 0 && bandBorders.second < targetLength);
    //printf("band: %d %d\n", bandBorders.first, bandBorders.second);

    Cell** matrix = new Cell*[targetLength];  // NOTE: First index is column, second is row.
    Cell* initialColumn = new Cell[queryLength];
    const int LOWER_SCORE_BOUND = INT_MIN + std::max(gapOpen, gapExt);
    for (int r = 0; r < queryLength; r++) {
        initialColumn[r].H = -1 * gapOpen - r * gapExt;
        initialColumn[r].E = LOWER_SCORE_BOUND;
    }

    Cell* prevColumn = initialColumn;
    int maxScore = INT_MIN;  // Max score so far, but only among cells that could be final.
    int H = INT_MIN;  // Current score.
    int c;
    for (c = 0; c < targetLength && maxScore < scoreLimit; c++) {
        matrix[c] = new Cell[queryLength];

        // First and last row in band for this column.
        int rBandStart = std::max(0, c - bandBorders.second);
        int rBandEnd = std::min(queryLength - 1, c + bandBorders.first);

        int uF, uH, ulH;
        if (rBandStart == 0) {
            uF = LOWER_SCORE_BOUND;
            uH = -1 * gapOpen - c * gapExt;
            ulH = c == 0 ? 0 : uH + gapExt;
        } else {
            uH = uF = LOWER_SCORE_BOUND;  // Out of band, so set to -inf.
            ulH = prevColumn[rBandStart - 1].H;
        }

        for (int r = rBandStart; r <= rBandEnd; r++) {
            int E = std::max(prevColumn[r].H - gapOpen, prevColumn[r].E - gapExt);
            int F = std::max(uH - gapOpen, uF - gapExt);
            int score = scoreMatrix[query[r] * alphabetLength + target[c]];
            H = std::max(E, std::max(F, ulH + score));
            /*
            printf("E: %d ", E);
            printf("F: %d ", F);
            printf("score: %d ", score);
            printf("ulH: %d ", ulH);
            printf("H: %d ", H);
            */

            // If mode is SW, track max score of all cells.
            // If mode is OV, track max score in last column.
            if (mode == OPAL_MODE_SW
                || (mode == OPAL_MODE_OV && c == targetLength - 1)) {
                maxScore = std::max(maxScore, H);
            }

            uF = F;
            uH = H;
            ulH = prevColumn[r].H;

            matrix[c][r].H = H;
            matrix[c][r].E = E;
            matrix[c][r].F = F;
        }

        // Set all cells that are out of band to -inf.
        for (int r = 0; r < rBandStart; r++) {
            matrix[c][r].E = matrix[c][r].H = matrix[c][r].F = LOWER_SCORE_BOUND;
        }
        for (int r = rBandEnd + 1; r < queryLength; r++) {
            matrix[c][r].E = matrix[c][r].H = matrix[c][r].F = LOWER_SCORE_BOUND;
        }

        if (mode == OPAL_MODE_HW || mode == OPAL_MODE_OV) {
            maxScore = std::max(maxScore, H);  // Track max score in last row.
        }
        prevColumn = matrix[c];
    }
    int lastColumnIdx = c - 1;

    result->startLocationTarget = 0;
    result->startLocationQuery = 0;
    result->scoreSet = 1;
    // Determine score and end location of alignment.
    switch (mode) {
    case OPAL_MODE_NW:
        opalSearchResultSetScore(result, H);
        result->endLocationTarget = targetLength - 1;
        result->endLocationQuery = queryLength - 1;
        break;
    case OPAL_MODE_HW:
        opalSearchResultSetScore(result, maxScore);
        result->endLocationTarget = lastColumnIdx;
        result->endLocationQuery = queryLength - 1;
        break;
    case OPAL_MODE_SW: case OPAL_MODE_OV:
        opalSearchResultSetScore(result, maxScore);
        result->endLocationTarget = lastColumnIdx;
        int r;
        for (r = 0; r < queryLength && matrix[lastColumnIdx][r].H != maxScore; r++);
        assert(r < queryLength);
        assert(matrix[lastColumnIdx][r].H == maxScore);
        result->endLocationQuery = r;
        break;
    default:
        assert(false);
    }

    // Construct alignment.
    // I reserve max size possibly needed for alignment.
    unsigned char* alignment = (unsigned char*) malloc(
        sizeof(unsigned char) * (result->endLocationQuery + result->endLocationTarget));
    int alignmentLength = 0;
    int rIdx = result->endLocationQuery;
    int cIdx = result->endLocationTarget;
    Cell::Field field = Cell::Field::H;  // Current field type.
    while (rIdx >= 0 && cIdx >= 0) {
        Cell cell = matrix[cIdx][rIdx];  // Current cell.

        // Determine to which cell and which field we should go next, move, and add operation to alignment.
        switch (field) {
        case Cell::Field::H:
            if (cell.H == cell.E) {
                field = Cell::Field::E;
            } else if (cell.H == cell.F) {
                field = Cell::Field::F;
            } else {
                alignment[alignmentLength++] = (query[rIdx] == target[cIdx] ? OPAL_ALIGN_MATCH
                                                : OPAL_ALIGN_MISMATCH);
                cIdx--; rIdx--;
            }
            break;
        case Cell::Field::E:
            field = (cell.E == matrix[cIdx - 1][rIdx].H - gapOpen) ? Cell::Field::H : Cell::Field::E;
            alignment[alignmentLength++] = OPAL_ALIGN_INS;
            cIdx--;
            break;
        case Cell::Field::F:
            field = (cell.F == matrix[cIdx][rIdx - 1].H - gapOpen) ? Cell::Field::H : Cell::Field::F;
            alignment[alignmentLength++] = OPAL_ALIGN_DEL;
            rIdx--;
            break;
        }
    }
    // I stop when matrix border is reached, so I have to add indels at start of alignment
    // manually (they do not have entry in operations). Only one of these two loops will trigger.
    while (rIdx >= 0) {
        alignment[alignmentLength] = OPAL_ALIGN_DEL;
        alignmentLength++; rIdx--;
    }
    while (cIdx >= 0) {
        alignment[alignmentLength] = OPAL_ALIGN_INS;
        alignmentLength++; cIdx--;
    }
    //printf("rIdx: %d, cIdx: %d\n", rIdx, cIdx);
    assert(rIdx == -1 && cIdx == -1);
    alignment = (unsigned char*) realloc(alignment, sizeof(unsigned char) * alignmentLength);
    revertArray(alignment, alignmentLength);
    // Store alignment to result.
    result->alignment = alignment;
    result->alignmentLength = alignmentLength;

    /*
    printf("Alignment: ");
    for (int j = 0; j < result->alignmentLength; j++)
        printf("%d ", result->alignment[j]);
    printf("\n");
    */

    // Cleanup
    delete[] initialColumn;
    for (int i = 0; i <= lastColumnIdx; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}



extern int opalSearchDatabase(
    unsigned char query[], int queryLength,
    unsigned char** db, int dbLength, int dbSeqLengths[],
    int gapOpen, int gapExt, int* scoreMatrix, int alphabetLength,
    OpalSearchResult* results[], const int searchType, int mode, int overflowMethod) {
// A.D.: this will not happen since SIMDe always takes care of correct SIMD extensions
#if 0
//#if !defined(__SSE4_1__) && !defined(__AVX2__)
//    return OPAL_ERR_NO_SIMD_SUPPORT;
#else
    // Calculate score and end location.
    int status;
    // Skip recalculation of already calculated sequences.
    bool *skip = new bool[dbLength];
    for (int i = 0; i < dbLength; i++) {
        skip[i] = (!opalSearchResultIsEmpty(*results[i])
                   && (searchType == OPAL_SEARCH_SCORE
                       || (results[i]->endLocationQuery >= 0 && results[i]->endLocationTarget >= 0)));
    }
    if (mode == OPAL_MODE_NW) {
        status = searchDatabase<OPAL_MODE_NW>(
            query, queryLength, db, dbLength, dbSeqLengths, gapOpen, gapExt,
            scoreMatrix, alphabetLength, results, searchType, skip, overflowMethod);
    } else if (mode == OPAL_MODE_HW) {
        status = searchDatabase<OPAL_MODE_HW>(
            query, queryLength, db, dbLength, dbSeqLengths, gapOpen, gapExt,
            scoreMatrix, alphabetLength, results, searchType, skip, overflowMethod);
    } else if (mode == OPAL_MODE_OV) {
        status = searchDatabase<OPAL_MODE_OV>(
            query, queryLength, db, dbLength, dbSeqLengths, gapOpen, gapExt,
            scoreMatrix, alphabetLength, results, searchType, skip, overflowMethod);
    } else if (mode == OPAL_MODE_SW) {
        status = searchDatabaseSW(
            query, queryLength, db, dbLength, dbSeqLengths,
            gapOpen, gapExt, scoreMatrix, alphabetLength,
            results, searchType, skip, overflowMethod);
    } else {
        status = OPAL_ERR_INVALID_MODE;
    }
    delete[] skip;
    if (status) return status;

    if (searchType == OPAL_SEARCH_ALIGNMENT) {
        // Calculate alignment of query with each database sequence.
        unsigned char* const rQuery = createReverseCopy(query, queryLength);
        for (int i = 0; i < dbLength; i++) {
            if (mode == OPAL_MODE_SW && results[i]->score == 0) {  // If it does not have alignment
                results[i]->alignment = NULL;
                results[i]->alignmentLength = 0;
                results[i]->startLocationQuery = results[i]->startLocationTarget = -1;
                results[i]->endLocationQuery = results[i]->endLocationTarget = -1;
            } else {
                //printf("%d %d\n", results[i]->endLocationQuery, results[i]->endLocationTarget);
                // Do alignment in reverse direction.
                int alignQueryLength = results[i]->endLocationQuery + 1;
                unsigned char* alignQuery = rQuery + queryLength - alignQueryLength;
                int alignTargetLength = results[i]->endLocationTarget + 1;
                unsigned char* alignTarget = createReverseCopy(db[i], alignTargetLength);
                OpalSearchResult result;
                findAlignment(
                    alignQuery, alignQueryLength, alignTarget, alignTargetLength,
                    gapOpen, gapExt, scoreMatrix, alphabetLength,
                    results[i]->score, &result, mode);
                //printf("%d %d\n", results[i]->score, result.score);
                assert(results[i]->score == result.score);
                // Translate results.
                results[i]->startLocationQuery = alignQueryLength - result.endLocationQuery - 1;
                results[i]->startLocationTarget = alignTargetLength - result.endLocationTarget - 1;
                results[i]->alignmentLength = result.alignmentLength;
                results[i]->alignment = result.alignment;
                revertArray(results[i]->alignment, results[i]->alignmentLength);
                free(alignTarget);
            }
        }
        free(rQuery);
    } else {
        for (int i = 0; i < dbLength; i++) {
            results[i]->alignment = NULL;
            results[i]->alignmentLength = -1;
            results[i]->startLocationQuery = -1;
            results[i]->startLocationTarget = -1;
        }
    }

    return 0;
#endif
}


extern int opalSearchDatabaseCharSW(
    unsigned char query[], int queryLength, unsigned char** db, int dbLength,
    int dbSeqLengths[], int gapOpen, int gapExt, int* scoreMatrix,
    int alphabetLength, OpalSearchResult* results[]) {
#if !defined(__SSE4_1__) && !defined(__AVX2__)
    return OPAL_ERR_NO_SIMD_SUPPORT;
#else
    bool* calculated = new bool[dbLength];
    for (int i = 0; i < dbLength; i++) {
        calculated[i] = false;
    }
    int resultCode = searchDatabaseSW_< SimdSW<char> >(
        query, queryLength, db, dbLength, dbSeqLengths, gapOpen, gapExt,
        scoreMatrix, alphabetLength, results, OPAL_SEARCH_SCORE, calculated,
        OPAL_OVERFLOW_SIMPLE);
    for (int i = 0; i < dbLength; i++) {
        if (!calculated[i]) {
            results[i]->score = -1;
            results[i]->scoreSet = 0;
        }
    }
    delete[] calculated;
    return resultCode;
#endif
}


extern void opalInitSearchResult(OpalSearchResult* result) {
    result->scoreSet = 0;
    result->startLocationTarget = result->startLocationQuery = -1;
    result->endLocationTarget = result->endLocationQuery = -1;
    result->alignment = NULL;
    result->alignmentLength = 0;
}

extern int opalSearchResultIsEmpty(const OpalSearchResult result) {
    return !result.scoreSet;
}

extern void opalSearchResultSetScore(OpalSearchResult* result, int score) {
    result->scoreSet = 1;
    result->score = score;
}
