/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

	prc.h: constants, structures, functions and macros used by PRC

                 $Revision: 1.47 $ $Date: 2005/09/16 16:33:27 $

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#ifndef __PRC_H__
#define __PRC_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

// release version
#define VERSION   \
"1.0.0"

// copyright notice printed out by all the binaries
#define COPYRIGHT \
"Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK"
#define COPYRIGHT2 \
"Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands"

// license printed out by all the binaries
#define LICENSE   \
"Freely distributed under the GNU General Public License"

// file formats known to all the binaries
#define KNOWN_EXTENSIONS                                         \
"  .prc     : binary PRC format\n"                               \
"  .mod     : ASCII or binary SAM3 model (one model per file)\n" \
"  .hmm     : ASCII HMMER2.0/3.0 model (one model per file)\n"   \
"  .fa      : single-sequence FASTA file (also .seq .fasta)\n"   \
"  .a2m     : multi-sequence aligned FASTA file (also .aln)\n"   \
"  .profile : binary PSI-BLAST profile (\"checkpoint file\")\n"  \
"             (also .prof .pbla .psi .psiblast)\n"

// maximum allowed profile HMM length
//
#define MAX_PROF_HMM_LENGTH 9999

// maximum and minimum doubles (used to determine overflows)
#define MIN_FLOAT 1.0e-35
#define MAX_FLOAT 1.0e35

// 'log(0)'
#define LOG_ZERO -9999.9

// the PRC log function^TM
#define LOG(x) ((x<MIN_FLOAT) ? LOG_ZERO : log(x))

// profile HMM topology:
//
// PLAN9  all 9 transitions between {M,I,D}, a la SAM (default)
// PLAN7  no I->D and D->I transitions, a la HMMER
//
// same for all object files in a build
//
#define PLAN7 0
#define PLAN9 1

#ifndef PROF_HMM_TRANS
#define PROF_HMM_TRANS PLAN9
#endif

// state space of the pair meta HMM:
//
// SPACE9  all states          =  MM, MI, IM, MD, DM, DD, ID, DI, II (default)
// SPACE8  all states bar II   =  MM, MI, IM, MD, DM, DD, ID, DI
// SPACE6  minimal space + DD  =  MM, MI, IM, MD, DM, DD
// SPACE5  minimal space       =  MM, MI, IM, MD, DM
//
// same for all object files in a build
//
#define SPACE9 0
#define SPACE8 1
#define SPACE6 2
#define SPACE5 3

#ifndef PAIR_HMM_STATES
#define PAIR_HMM_STATES SPACE9
#endif

// transitions in the pair meta HMM:
//
// ALL_TRANS  all transitions allowed by the profile HMM (default)
// VIA_MM     transitions must start/end in the MM state
//
// same for all object files in a build
//
#define ALL_TRANS 0
#define VIA_MM    1

#ifndef PAIR_HMM_TRANS
#define PAIR_HMM_TRANS ALL_TRANS
#endif

// implementation of the search algorithms:
//
// LINSPACE  mutliply odds (default)
// LOGSPACE  add log odds
//
// DP_search.c is present in both versions in a build; the LINSPACE version
// calls the LOGSPACE version whenever it is about to overflow (see DP_search.c
// for details)
//
#define LINSPACE 0
#define LOGSPACE 1

#ifndef ALGO_SPACE
#define ALGO_SPACE LINSPACE
#endif

#if     ALGO_SPACE == LOGSPACE
#define LIN_LOG(x) LOG(x)
#define ZERO       LOG_ZERO
#define ONE        0.0
#define x          +
#define div        -
#define xeq        +=
#elif   ALGO_SPACE == LINSPACE
#define LIN_LOG(x) (x)
#define ZERO       0.0
#define ONE        1.0
#define x          *
#define div        /
#define xeq        *=
#endif

/*
// more implementation settings:
//
// NORM_TRANS==0  don't normalize transition Ps in the pair HMM to 1
// NORM_TRANS==1  normalize them
//
// as far as I can see the normalization is likely to be costly, and there's
// unlikely to be much benefit from it, but .... I should try it one day
//
#define NORM_TRANS 0
*/

// specific dynamic programming algorithms:
//
// VITERBI    standard Viterbi
// VIT_BACK   is to Viterbi as backward is to forward
// FORWARD    standard forward algorithm
// BACKWARD   standard backward algorithm
// OPT_ACCUR  the Holmes/Durbin optimal accuracy algorithm
//
#define  VITERBI      0
#define  VIT_BACK     1
#define  FORWARD      2
#define  BACKWARD     3
#define  OPT_ACCUR    4

#ifndef ALGORITHM
#define ALGORITHM FORWARD
#endif

// what to do about paths:
//
// SUM_PATHS  sum across all paths (forward/backward)
// BEST_PATH  best path (Viterbi/Viterbi_back)
//
// DP.c is present in both versions in a build; which gets called at runtime
// depends on program parameters
//
#define SUM_PATHS 0
#define BEST_PATH 1

// search direction:
//
// DIR_FORWARD   forward direction
// DIR_BACKWARD  backward direction
//
#define DIR_FORWARD  0
#define DIR_BACKWARD 1

#if     ALGORITHM == VITERBI
#define DIRECTION   DIR_FORWARD
#define ALGO_PATHS  BEST_PATH
#elif   ALGORITHM == VIT_BACK
#define DIRECTION   DIR_BACKWARD
#define ALGO_PATHS  BEST_PATH
#elif   ALGORITHM == FORWARD
#define DIRECTION   DIR_FORWARD
#define ALGO_PATHS  SUM_PATHS
#elif   ALGORITHM == BACKWARD
#define DIRECTION   DIR_BACKWARD
#define ALGO_PATHS  SUM_PATHS
#elif   ALGORITHM == OPT_ACCUR
#define DIRECTION   DIR_FORWARD
#define ALGO_PATHS  BEST_PATH
#endif

// program parameters:
//
// PARAM_FORW_BACK  run forward/backward + Holmes/Durbin
// PARAM_VITERBI    run Viterbi + Viterbi traceback
//
#define PARAM_FORW_BACK 0
#define PARAM_VITERBI   1

// scoring functions
#define  DOT1   0
#define  DOT2   1

// alignment modes
//
#define  LOCAL  0
#define  GLOBAL 1

// alignment styles
#define  ALIGN_NONE  0
#define  ALIGN_PRC1  1
#define  ALIGN_PRC2  2
#define  ALIGN_SAM1  3
#define  ALIGN_SAM2  4


// HMM information necessary for printing out search results
//
// (in memory for the whole library, so this should be kept to a minimum)
//
struct HMMinfo {
	char   *filename;        // which file did the profile HMM come from?
	char   *name;            // name of the profile HMM
	int    M;                // number of non-zero match states ("model length")
};
typedef struct HMMinfo HMMinfo;

// the profile HMM under comparison
//
// PRC follows the SAM convention that Pt[i][tXY] is the probability of
// transition *to* state Y_i, from either X_i or X_i-1
//
// segments 0 and M+3 are padding and everything connected with them is set to
// zero; M_1 is Begin, M_(M+2) is End, so M_2 is the first match state, etc.
//
// so, when printing out regions for the user, the segment number must be
// decreased by one (because our M_m is actually match m-1 in the HMM)
//
struct HMM {
	HMMinfo *i;              // see above

	double  **Pmat;           // Pmat[2..M+1][0..19],   match emission Ps
	double  **Pins;           // Pins[1..M+1][0..19],   insert emission Ps
	double  **Pt;             //   Pt[1..M+2][0..8|10], transition Ps

	double  **Smat;           // Pmat[i]/sqrt(hmm1->PinsJ[i] * hmm2->PinsJ[i])
	double  **Sins;           // Pins[i]/sqrt(hmm1->PinsJ[i] * hmm2->PinsJ[i])
	double  **St_lin;         // transitions scores; variously normalized wrt PtJJ
	double  **St_log;         // St_log = log(St_lin);

	double  PinsJ[20];        // PinsJ[0..19], simple null emission Ps (state J)
	double  PtJJ;             // J->J transition P

	unsigned char *discrete_alphabet; //the match states transalted into a discrete alphabet
};
typedef struct HMM HMM;

typedef struct {
	int size;        // The number of context profiles stored in the library
	int wlen;        // The window size either 13 for cs4000 of 1 for the discrete alphabet cs219
	int nalph;       // The alphabet size, currently 20 AA long
	double pc_admix;
	double pc_ali;
	double weight_center;
	double weight_decay;
	double weight_as;
	int center;
	double *weights; // The weights for calculating the posterior probabilities (the middle of the sequence is considered more important)
	double *prior;   // The prior probabilities of the context profiles 
	double ***probs; // The probabilities of the amino acids stored in a matrix of dimension: [size][wlen][nalp]
} CONTEXT_LIB;

#if     ALGO_SPACE == LOGSPACE
#define St St_log
#elif   ALGO_SPACE == LINSPACE
#define St St_lin
#endif

// indices for state transitions, stored in Pt[][] & St_xxx[][]
//
#define  tBM  0
#define  tME  1
#define  tMM  2
#define  tMI  3
#define  tMD  4
#define  tIM  5
#define  tII  6
#define  tDM  7
#define  tDD  8

#if      PROF_HMM_TRANS==PLAN9
#define  tID  9
#define  tDI  10
#endif

#if     PROF_HMM_TRANS == PLAN7
#define N_PROF_HMM_TRANS   9
#elif   PROF_HMM_TRANS == PLAN9
#define N_PROF_HMM_TRANS   11
#endif


// region in a DP matrix, understood as
//
//    startX <= mX <= endX
//
struct REGION {
	int    start1, end1, start2, end2;
};
typedef struct REGION REGION;

// match between two HMMs
//
struct MATCH {
	HMMinfo *hmm1i, *hmm2i;  // infos of the two HMMs that were compared
	int    free_hmm2i;       // free hmm2i when deallocating this match?
	int    match_n;          // match number (when more than 1 between two HMMs)
	REGION proper;           // the matching region (first and last proper state)
	REGION in_out;           // the matching region (InIn and OutOut)

	double  sum;              // only for FORWARD: the overall score
	double  simple;           // the score due to the match region
	double  reverse;          // same as above, minus the reverse score
	double E_value;          // the E-value calculated from the reverse score

	int    n_trace;          // length of the match region
	int    **trace;          // trace[0..n_trace-1][tr_XX], see below for tr_XX
	char   *alignment;       // the full alignment string

	struct MATCH *next;      // next match, if any
};
typedef struct MATCH MATCH;

// indices for **trace
#define  tr_m1     0
#define  tr_m2     1
#define  tr_state  2


// program parameters and shared values
//
struct PARAMS {
	int    algorithm;        // PARAM_FORW_BACK or PARAM_VITERBI
	int    MM_function;      // DOT1 or DOT2
	int    align_mode1;      // LOCAL or GLOBAL = alignment modes to HMM1
	int    align_mode2;      // same for HMM2
	int    align_style;      // ALIGN_NONE, ALIGN_PRC1, ALIGN_PRC2, ALIGN_SAM1, ALIGN_SAM2
	int    E_values;         // are we going to calculate E-values?
	double Emax;             // only report hits with E-value <= Emax
	double stop;             // stop looking for more hits when simple < stop
	long   int max_hits;     // stop looking for more hits when match > max_hits

	char   *output_scores;   // file to which save the main output
	char   *output_aligns;   // file to which save alignments
	char   *model_file1;     // file from which to read the query
	char   *model_file2;     // file from which to read the target(s)
	int    library;          // is the 2nd argument a library (or a model)?
	int    tree;             // should we build a distance matrix and a UPMGA tree?
	int    filt;             // should we filter the library using a fast discrete alpabet smith-waterman?
	char   *output_distmat;  // file to which save the distance matrix
	char   *output_tree;     // file to which save the UPGMA tree
	long   int n_models;     // number of models in the library (if library)

	double p[3];             // parameters of the E-value distribution
	MATCH  **matches;        // matches to be printed out
	long   int n_matches;    // size of the ^^^ array
	long   int n_unrel;      // number of unrelated matches, used for E-values
};
typedef struct PARAMS PARAMS;

// is m a weak match?
#define WEAK_MATCH(p,m) \
( (m->simple >= p->stop) && !((p->E_values) && (m->E_value > 10.0*p->Emax)) )

// is m a significant match?
#define SIG_MATCH(p,m) \
( (m->simple >= p->stop) && !((p->E_values) && (m->E_value > p->Emax)) )

// data used by the core dynamic programming algorithms
//
struct DPdata {
	HMM    *hmm1, *hmm2;     // the two models under comparison
	PARAMS *p;               // program parameters
	REGION region;           // which region of the matrix to search
	int    algorithm;        // FORWARD or VITERBI or BACKWARD or ...

	int    X_status;         // X_NOT_SET, X_SET or X_LIN2LOG
	double  ***X;             // [1..M1+1][1..M2+1][0..2|3] = M-M, M-I etc. scores
	double  ***S1;            // [0..M1+2][0..M2+2][0..4|5|7] = the DP matrix
	double  ***S2;            // S2 is like S1 but is used for BACKWARD
	int    ***S_tr;          // matrix for logging the best choice (VIT, OPT_ACC)

	int    lin_overflow;     // was there an over/underflow in LINSPACE?
	double  sum_S;            // sum of vvv (rather than the best)
	double  best_S;           // best S[sInIn] x StMMEE (or S[sOutOut] x StBBMM)
	int    best_m1, best_m2; //  & where it occurred

	MATCH  *list;            // the last match; the next one is list->next etc.
};
typedef struct DPdata DPdata;

// status of the X matrix
#define  X_NOT_SET   0
#define  X_SET       1
#define  X_LIN2LOG   2

// indices for the X matrix (match-match, match-insert etc. scores)
//
#if      PAIR_HMM_STATES != SPACE9
#define  xMM       0
#define  xMI       1
#define  xIM       2
#elif    PAIR_HMM_STATES == SPACE9
#define  xMM       0
#define  xMI       1
#define  xIM       2
#define  xII       3
#endif

// number of indices above
//
#if      PAIR_HMM_STATES != SPACE9
#define  N_X_STATES  3
#elif    PAIR_HMM_STATES == SPACE9
#define  N_X_STATES  4
#endif

// indices for the S1/S2 matrix:
//
#if      PAIR_HMM_STATES == SPACE5
#define  sMM       0
#define  sMI       1
#define  sIM       2
#define  sMD       3
#define  sDM       4
#define  sBE       5
#elif    PAIR_HMM_STATES == SPACE6
#define  sMM       0
#define  sMI       1
#define  sIM       2
#define  sMD       3
#define  sDM       4
#define  sDD       5
#define  sBE       6
#elif    PAIR_HMM_STATES == SPACE8
#define  sMM       0
#define  sMI       1
#define  sIM       2
#define  sMD       3
#define  sDM       4
#define  sDD       5
#define  sID       6
#define  sDI       7
#define  sBE       8
#elif    PAIR_HMM_STATES == SPACE9
#define  sMM       0
#define  sMI       1
#define  sIM       2
#define  sII       3
#define  sMD       4
#define  sDM       5
#define  sDD       6
#define  sID       7
#define  sDI       8
#define  sBE       9
#endif

// N.B. sBE isn't actually stored in the DP matrix
//
#if     PAIR_HMM_STATES == SPACE5
#define N_PAIR_HMM_STATES   5
#elif   PAIR_HMM_STATES == SPACE6
#define N_PAIR_HMM_STATES   6
#elif   PAIR_HMM_STATES == SPACE8
#define N_PAIR_HMM_STATES   8
#elif   PAIR_HMM_STATES == SPACE9
#define N_PAIR_HMM_STATES   9
#endif


// SORT THIS BIT OUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// function definitions
// ... from HMM.c
HMMinfo* free_HMM(HMM*);
void free_HMMinfo(HMMinfo*);
void set_HMM_align_mode(HMM*, int, int, int, int);
void scorify_HMM_M_I(HMM*, double*);
void print_HMM(FILE*, HMM*);
HMM* reverse_HMM(HMM*);
void write_HMM_PRC_binary(FILE*, HMM*);
void write_HMM_PRC_binary_with_discrete_alphabet(FILE*, HMM*);
HMM* read_HMM(char*);
int known_extension(char*);

// ... from DP_search.c
void search_LINSPACE_VITERBI(DPdata*, REGION*);
void search_LOGSPACE_VITERBI(DPdata*, REGION*);
void search_LINSPACE_VIT_BACK(DPdata*, REGION*);
void search_LOGSPACE_VIT_BACK(DPdata*, REGION*);
void search_LINSPACE_FORWARD(DPdata*, REGION*);
void search_LOGSPACE_FORWARD(DPdata*, REGION*);
void search_LINSPACE_BACKWARD(DPdata*, REGION*);
void search_LOGSPACE_BACKWARD(DPdata*, REGION*);
void search_LINSPACE_OPT_ACCUR(DPdata*, REGION*);
void search_LOGSPACE_OPT_ACCUR(DPdata*, REGION*);

// ... from DP_traceback.c
void traceback_VITERBI(DPdata*, int, int);
void traceback_OPT_ACCUR(DPdata*, int, int);

// ... from DP_common.c
void malloc_DPdata_X(DPdata*);
void malloc_DPdata_S2(DPdata*);
void malloc_DPdata_S_tr(DPdata*);
MATCH* free_DPdata(DPdata*);
void write_Sxx_gnuplot(char*, DPdata*, double***, int);
void print_DP(DPdata*, double***, int);
MATCH* run_two_HMMs(PARAMS*, HMM*, HMM*, HMM*, int, int);
MATCH* simple_run_two_HMMs(PARAMS*, HMM*, HMM*, int, int);

// ... from match.c
MATCH* free_match(MATCH*);
void free_matches(MATCH*);
int count_matches(MATCH*);
MATCH* join_lists(MATCH*, MATCH*);
int copy_sort_matches(PARAMS*, MATCH*);

// ... from region.c
void copy_region(REGION*, REGION*);
void set_region(REGION*, int, int, int, int);

// .. from stats.c
double E_value(long int, double, double*);
char* nice_string(double);
void fit(PARAMS*, MATCH*);

// .. from context.c
CONTEXT_LIB *read_context_library(char *filename);
double **sequence_to_profile(char *sequence, CONTEXT_LIB *cs4000);
unsigned char *profile_to_discrete_alphabet(double **profile, int profile_length, CONTEXT_LIB *cs219);
void free_context_library(CONTEXT_LIB *context_lib);
int *assign_match_columns_by_gap_rule(char **alignment, int number_of_sequences, double gap_threshold);
double *position_specific_weights_and_diversity(char **alignment, int number_of_sequences, double ***pdwm_in);
double *del_position_specific_weights_and_diversity(char **alignment, int number_of_sequences, double ***pdwm_in);
double *ins_position_specific_weights_and_diversity(char **alignment, int number_of_sequences, double ***pdwm_in);
double global_weights_and_diversity(char **alignment, int number_of_sequences, double **global_weights_in, bool neff_sum_pairs);
char **extract_match_columns(char **alignment, int nr_sequences, int *match_column);
char **extract_insert_columns(char **alignment, int nr_sequences, int *match_column);
char **extract_deletion_columns(char **alignment, int nr_sequences, int *match_column);
double **alignment_to_count_profile(char **alignment, int number_of_sequences, double **neff, bool pos_weight, bool neff_sum_pairs);
double **count_profile_to_profile(double **count_profile, int profile_length, double *neff, CONTEXT_LIB *cs4000);
int make_profile_transitions_and_neff(char **alignment, int nr_sequences, double ***profile_in, double ***transitions_in, double ***NeffMID_in);

// the dieX() series of macros
#define  die() { fprintf(stderr, \
                 "\n\nERROR ON LINE %d OF FILE %s, EXITING!\n\n", \
	         __LINE__, __FILE__); fflush (stderr); exit(1); }
#define  die1(a)         { fprintf(stderr, "\n"); \
                           fprintf(stderr, a);             die(); }
#define  die2(a,b)       { fprintf(stderr, "\n"); \
                           fprintf(stderr, a, b);          die(); }
#define  die3(a,b,c)     { fprintf(stderr, "\n"); \
                           fprintf(stderr, a, b, c);       die(); }
#define  die4(a,b,c,d)   { fprintf(stderr, "\n"); \
                           fprintf(stderr, a, b, c, d);    die(); }
#define  die5(a,b,c,d,e) { fprintf(stderr, "\n"); \
                           fprintf(stderr, a, b, c, d, e); die(); }

// file I/O
#define  open_file_or_die(stream, name, how) \
         { if( (stream = fopen( name, how )) == NULL ) \
           { char error[500]; \
             sprintf(error, "\nError opening file '%s'", name); \
	     perror(error); die(); }; }
#define  fgets_null_check(fgets, filename) \
         { if( (fgets) == NULL ) \
           { fprintf(stderr, "\nError reading '%s': unexpected end of file\n",\
                     filename); \
             die(); }; }

// memory allocation macros
#define  malloc_1D_array(a,t,x) \
         { if((a=(t*)malloc(sizeof(t)*(x)))==NULL) \
           die2("\nmalloc of %ldB failed\n", sizeof(t)*((long int)(x))); }
//           printf("# malloc_1D_array on line %d of file %s: %lx from %p\n",
//		  __LINE__, __FILE__, sizeof(t)*((long int)(x)), a); }
#define  calloc_1D_array(a,t,x) \
         { if((a=(t*)calloc(x,sizeof(t)))==NULL) \
           die2("\ncalloc of %ldB failed\n", sizeof(t)*((long int)(x))); }
//           printf("# calloc_1D_array on line %d of file %s: %lx from %p\n",
//		  __LINE__, __FILE__, sizeof(t)*((long int)(x)), a); }
#define  malloc_2D_array(a,t,x,y) \
         { int _i; malloc_1D_array(a,t*,x); malloc_1D_array(a[0],t,(x)*(y)); \
           for(_i=1; _i<(x); _i++) a[_i]=a[0]+_i*(y); }
#define  calloc_2D_array(a,t,x,y) \
         { int _i; calloc_1D_array(a,t*,x); calloc_1D_array(a[0],t,(x)*(y)); \
           for(_i=1; _i<(x); _i++) a[_i]=a[0]+_i*(y); }
#define  malloc_3D_array(a,t,x,y,z) \
         { int _x, _y; malloc_2D_array(a,t*,x,y); \
           malloc_1D_array(a[0][0],t,(x)*(y)*(z)); \
           for(_x=0;_x<(x);_x++) for(_y=0;_y<(y);_y++) \
           a[_x][_y] = a[0][0] + (z)*(_y + (y)*_x); }
#define  calloc_3D_array(a,t,x,y,z) \
         { int _x, _y; malloc_2D_array(a,t*,x,y); \
           calloc_1D_array(a[0][0],t,(x)*(y)*(z)); \
           for(_x=0;_x<(x);_x++) for(_y=0;_y<(y);_y++) \
           a[_x][_y] = a[0][0] + (z)*(_y + (y)*_x); }
#define  free_unless_null(a) \
         { if(a!=NULL){ free(a); a=NULL; } }
#define  free_2D_array(a) \
         { if(a!=NULL){ free(a[0]); free(a); a=NULL; } }
#define  free_3D_array(a) \
         { if(a!=NULL){ free(a[0][0]); free(a[0]); free(a); a=NULL; } }

#endif /* #ifndef __PRC_H__ */
