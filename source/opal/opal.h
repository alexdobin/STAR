#ifndef OPAL_H
#define OPAL_H

/***********************************************************************************
 *  - adapts score type and reports overflow if score does not fit in int
 *  - calculating column by column (not 4 of them at once)
 *  - db sequences are not padded
 *  - using saturation arithmetic when possible
 *  - works for SSE4.1 and higher
 *************************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

// Error codes
#define OPAL_ERR_OVERFLOW 1 //!< Returned when score overflow happens. Happens only if score can't fit in int.
#define OPAL_ERR_NO_SIMD_SUPPORT 2 //!< Returned if available SIMD is not SSE4.1 or higher.
#define OPAL_ERR_INVALID_MODE 3 //!< Returned when given mode is invalid.

// Modes
#define OPAL_MODE_NW 0
#define OPAL_MODE_HW 1
#define OPAL_MODE_OV 2
#define OPAL_MODE_SW 3

// Overflow handling
#define OPAL_OVERFLOW_SIMPLE 0
#define OPAL_OVERFLOW_BUCKETS 1

// Search types
#define OPAL_SEARCH_SCORE 0 //!< Search finds only score -> it is the fastest search.
#define OPAL_SEARCH_SCORE_END 1 //!< Search finds score and end location of alignment.
#define OPAL_SEARCH_ALIGNMENT 2 //!< Search finds score, start and end location of alignment and alignment.

// Alignment operations
#define OPAL_ALIGN_MATCH 0  //!< Match.
#define OPAL_ALIGN_DEL 1  //!< Deletion from query (insertion to target).
#define OPAL_ALIGN_INS 2  //!< Insertion to query (deletion from target).
#define OPAL_ALIGN_MISMATCH 3  //!< Mismatch.

    /**
     * Contains score and alignment information.
     * If there are multiple possible alignments, one whose end
     * has smallest position in target and then smallest position in row is used.
     */
    struct OpalSearchResult {
        int scoreSet;  //!< If 1, result contains at least score. If 0, whole result is empty.
        int score;

        //!< 0-indexed position in target where aligment ends. -1 if not set.
        int endLocationTarget;
        //!< 0-indexed position in query where aligment ends. -1 if not set.
        int endLocationQuery;

        //!< 0-indexed position in target where alignment starts. -1 if not set.
        int startLocationTarget;
        //!< 0-indexed position in query where alignment starts. -1 if not set.
        int startLocationQuery;

        /**
         * Alignment is sequence of operations:
         *  - OPAL_ALIGN_MATCH stands for match.
         *  - OPAL_ALIGN_DEL stands for deletion from query (insertion to target).
         *  - OPAL_ALIGN_INS stands for insertion to query (deletion from target).
         *  - OPAL_ALIGN_MISMATCH stands for mismatch.
         * Alignment aligns query to target from begining of query till end of query.
         * If gaps are not penalized, they are not in alignment.
         * Needed memory is allocated and given pointer is set to it.
         * Important: Do not forget to free memory allocated for alignment! Use free().
         */
        unsigned char* alignment; //!< NULL if there is no alignment.
        int alignmentLength; //!< 0 if there is no alignment.
    };

    /**
     * Initializes result to empty result.
     * Use this to initialize result object for first time, or to reset it.
     * If reseting it, make sure you free alignment if you do not need it anymore.
     * @param result
     */
    void opalInitSearchResult(OpalSearchResult* result);

    /**
     * @return 1 if result is empty, 0 if not (has at least score).
     */
    int opalSearchResultIsEmpty(const OpalSearchResult result);

    void opalSearchResultSetScore(OpalSearchResult* result, int score);


    /**
     * Compares query sequence with each database sequence and returns similarity scores.
     * Uses one of following alignment algorithms in combination with afine gaps(GOTOH):
     *   SW, NW, HW, OV.
     * Sequences are not represented as arrays of letters, but as arrays of indices of
     * letters in alphabet. For example, if alphabet is {A,C,T,G} and sequence is ACCTCAG
     * it will be represented as 0112103.
     * Opening of gap is penalized with gapOpen, while gap extension is penalized with
     * gapExt. Therefore, gap of length n will have penalty of gapOpen + (n - 1) * gapExt.
     * gapOpen, gapExt and scores from scoreMatrix must be in (INT_MIN/2, INT_MAX/2).
     * Detects overflow only for SW!
     * Although not crucial, sorting database sequences ascending by length can give
     * some extra speed.
     * @param [in] query Query sequence.
     * @param [in] queryLength Length of query sequence.
     * @param [in] db Array of database sequences (each sequence is also an array).
     * @param [in] dbLength Number of database sequences.
     * @param [in] dbSeqLengths Array of lengths of database sequences.
     * @param [in] gapOpen  Non-negative penalty for gap opening.
     * @param [in] gapExt  Non-negative penalty for gap extension.
     * @param [in] scoreMatrix Matrix of dimensions (alphabetLength, alphabetLength).
     *     It is array of length alphabetLength * alphabetLength, where memory is organized
     *     row by row: row0row1row2...rowN.
     *     When there is a (mis)match of element Q from query and element T from target,
     *     its score is read from scoreMatrix[Q * alphabetLength + T].
     * @param [in] alphabetLength
     * @param [in|out] results Results of search are written here, for each sequence.
     *     If a result has score and end location already calculated, they will not be
     *     calculated again, and if you are searching for alignment, they will be used
     *     to find alignment (that way you can reuse previous result).
     *     If you provide score and end location that are incorrect, behavior is undefined.
     * @param [in] searchType Defines what type of search will be conducted.
     *     OPAL_SEARCH_SCORE: find score.
     *     OPAL_SEARCH_SCORE_END: find score and end location.
     *     OPAL_SEARCH_ALIGNMENT: find score, end location, start location and alignment.
     *       Finding of alignment takes significant amount of time and memory.
     * @param [in] mode Mode of alignment, different mode means different algorithm.
     *     OPAL_MODE_NW: global alignment (Needleman-Wunsch)
     *     OPAL_MODE_HW: semi-global. Gap at query start and gap at query end
     *                    are not penalized.
     *                                 DBSEQ
     *                                _QUERY_
     *     OPAL_MODE_OV: semi-global. Gap at query start, gap at query end,
     *                    gap at dbseq start and gap at dbseq end are not penalized.
     *                                _DBSEQ_
     *                                _QUERY_
     *     OPAL_MODE_SW: local alignment (Smith-Waterman)
     * @param [in] overflowMethod Method that defines behavior regarding overflows.
     *     OPAL_OVERFLOW_SIMPLE: all sequences are first calculated using
     *         char precision, those that overflowed are then calculated using
     *         short precision, and those that overflowed again are calculated
     *         using integer precision.
     *     OPAL_OVERFLOW_BUCKETS: database is divided into buckets and each
     *         bucket is calculated independently. When overflow occurs,
     *         calculation is resumed with higher precision for all
     *         following sequences in that bucket.
     * @return 0 if all okay, error code otherwise.
     */
    int opalSearchDatabase(
        unsigned char query[], int queryLength, unsigned char* db[], int dbLength,
        int dbSeqLengths[], int gapOpen, int gapExt, int* scoreMatrix,
        int alphabetLength, OpalSearchResult* results[],
        const int searchType, int mode, int overflowMethod);

    /**
     * Same like opalSearchDatabase, with few small differences:
     * - uses char for score representation
     * - works in OPAL_OVERFLOW_SIMPLE mode: does not stop on overflow
     * - if sequence i overflows, sequence[i] is set to -1
     */
    int opalSearchDatabaseCharSW(
        unsigned char query[], int queryLength, unsigned char** db, int dbLength,
        int dbSeqLengths[], int gapOpen, int gapExt, int* scoreMatrix,
        int alphabetLength, OpalSearchResult* results[]);

#ifdef __cplusplus
}
#endif

#endif /* OPAL_H */
