/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

	      DP_search.c: the core dynamic programming routines

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prc.h"


/*
  This file gets "alternatively spliced" to produce all the various
  DP_search_xxx.o files


  === WARNING: THIS IS A BIT OUT OF DATE ===

  An overview of the LINSPACE vs LOGSPACE issue
  - - - - - - - - - - - - - - - - - - - - - - -

  The dymamic programming algorithms in this file spend most of their time
  dealing with doubles in linear space. This makes the implementation simple &
  fast (adding two probabilities becomes trivial), but also susceptible to
  overflows. Underflows are less of a problem, because I am dealing with ratios
  (model P divided by null P) rather than raw probabilities, so an underflow
  automatically implies a poorly matching region ... and I'm not particularly
  fussed about their exact scores.

  Since empirically the overflows only happen for very closely related profiles
  (usually ones where the homology is obvious using even simple pairwise
  methods), and because such pairs are rare, when I encounter one I can afford
  to branch the code and spend a comparatively long period of time dealing with
  it seperately.

  :

*/

// calculate the score between two emission states

#if ALGORITHM != OPT_ACCUR
static double score_dot1(double *Pemm1, double *Pemm2, double null_score)
{
	double sum = 0.0;
	int   i;

	for (i = 0; i < 20; i++) {
		sum += Pemm1[i] * Pemm2[i];
	}

	return LIN_LOG(sum) div null_score;
}

static double score_dot2(double *Semm1, double *Semm2)
{
	double sum = 0.0;
	int   i;

	for (i = 0; i < 20; i++) {
		sum += Semm1[i] * Semm2[i];
	}

	return LIN_LOG(sum);
}
#endif


//
//
#if   ALGO_SPACE==LINSPACE

#if   ALGORITHM == VITERBI
#define search search_LINSPACE_VITERBI
#elif ALGORITHM == VIT_BACK
#define search search_LINSPACE_VIT_BACK
#elif ALGORITHM == FORWARD
#define search search_LINSPACE_FORWARD
#elif ALGORITHM == BACKWARD
#define search search_LINSPACE_BACKWARD
#elif ALGORITHM == OPT_ACCUR
#define search search_LINSPACE_OPT_ACCUR
#endif

#elif ALGO_SPACE==LOGSPACE

#if   ALGORITHM == VITERBI
#define search search_LOGSPACE_VITERBI
#elif ALGORITHM == VIT_BACK
#define search search_LOGSPACE_VIT_BACK
#elif ALGORITHM == FORWARD
#define search search_LOGSPACE_FORWARD
#elif ALGORITHM == BACKWARD
#define search search_LOGSPACE_BACKWARD
#elif ALGORITHM == OPT_ACCUR
#define search search_LOGSPACE_OPT_ACCUR
#endif

#endif


// comparison functions to make the search() code more tractable

#include "auto_DP_which_comps.c"
#include "auto_DP_comps.c"


//#define DEBUG_SEARCH
//#define DEBUG_BOUNDARY
void search(DPdata* dp, REGION *r)
{
	HMM   *hmm1 = dp->hmm1, *hmm2 = dp->hmm2;
	int   m1_start = r->start1, m1_end = r->end1;
	int   m2_start = r->start2, m2_end = r->end2;
	double **St1 = hmm1->St, **St2 = hmm2->St, ***S;
	double StBBMM, StMMEE, temp;

#if ALGORITHM != OPT_ACCUR
	double null_score = score_dot1(hmm1->PinsJ, hmm2->PinsJ, ONE);
	double ***X;
#endif

#if (ALGORITHM == FORWARD) || (ALGORITHM == BACKWARD)
	double ***other_S;
#endif

#if ALGO_PATHS == BEST_PATH
	int   ***Tr;
#endif

	int   m1, m2, i;

#ifdef DEBUG_SEARCH
	char  *algo[]  = { "VITERBI", "VIT_BACK",
	                   "FORWARD", "BACKWARD", "OPT_ACCUR"
	                 };
	char  *space[] = { "LINSPACE", "LOGSPACE" };
	char  *paths[] = { "SUM_PATHS", "BEST_PATH" };
	char  *dir[]   = { "DIR_FORWARD", "DIR_BACKWARD" };

	printf("\n\nSEARCH:\n\n");
	printf("  hmm1: %s\n", hmm1->i->name);
	printf("  hmm2: %s\n", hmm2->i->name);
	printf("  algo: %s\n", algo[ALGORITHM]);
	printf(" space: %s\n", space[ALGO_SPACE]);
	printf(" paths: %s\n", paths[ALGO_PATHS]);
	printf("   dir: %s\n", dir[DIRECTION]);
	printf("  reg1: %d - %d\n", m1_start, m1_end);
	printf("  reg2: %d - %d\n", m2_start, m2_end);
	printf("\n");
#endif

#if ALGORITHM != OPT_ACCUR
	if ((dp->X_status == X_NOT_SET) && (dp->X == NULL)) {
		malloc_DPdata_X(dp);
	}
	X = dp->X;
#endif

#if ALGO_PATHS == SUM_PATHS
	malloc_DPdata_S2(dp);
#endif

#if ALGO_PATHS == BEST_PATH
	malloc_DPdata_S_tr(dp);
	Tr = dp->S_tr;
#endif

#if   ALGORITHM == VITERBI
	S = dp->S1;
#elif ALGORITHM == VIT_BACK
	S = dp->S2;
#elif ALGORITHM == FORWARD
	S = dp->S1;
	other_S = dp->S2;
#elif ALGORITHM == BACKWARD
	S = dp->S2;
	other_S = dp->S1;
#elif ALGORITHM == OPT_ACCUR
	// N.B. dp->algorithm is the last algorithm that was run
	if (dp->algorithm == FORWARD) {
		S = dp->S2;
	} else if (dp->algorithm == BACKWARD) {
		S = dp->S1;
	} else {
		die4("Can only run OPT_ACCUR after FORWARD=#%d or BACKWARD=#%d,\n"
		     "not ALGORITHM=#%d, whatever it is!",
		     FORWARD, BACKWARD, ALGORITHM);
	};
#endif

	dp->algorithm = ALGORITHM;

	// don't touch best_xx with OPT_ACCUR!
#if ALGORITHM != OPT_ACCUR
	dp->best_S  = ZERO;
	dp->best_m1 = 0;
	dp->best_m2 = 0;
#endif
#if ALGO_PATHS == SUM_PATHS
	dp->sum_S = ZERO;
#endif


#if   DIRECTION == DIR_FORWARD
#define StIn         StBBMM
#define StOut        StMMEE
#define m1_00_row    m1_start-1
#define m2_00_col    m2_start-1
#define m1_11_row    m1_start
#define m2_11_col    m2_start
#define m1_diff      m1++
#define m2_diff      m2++
#define m1_end_cond  m1<=m1_end
#define m2_end_cond  m2<=m2_end
#elif DIRECTION == DIR_BACKWARD
#define StIn         StMMEE
#define StOut        StBBMM
#define m1_00_row    m1_end+1
#define m2_00_col    m2_end+1
#define m1_11_row    m1_end
#define m2_11_col    m2_end
#define m1_diff      m1--
#define m2_diff      m2--
#define m1_end_cond  m1>=m1_start
#define m2_end_cond  m2>=m2_start
#endif

	// set up the boundary
#ifdef DEBUG_BOUNDARY
	printf("setting up the boundary:\n\n");
#endif

	for (m1 = m1_00_row; m1_end_cond; m1_diff) {
#ifdef DEBUG_BOUNDARY
		printf("S[%d][%d][sXX]=ZERO\n", m1, m2_00_col);
#endif

		for (i = 0; i < N_PAIR_HMM_STATES; i++) {
			S[m1][m2_00_col][i] = ZERO;
		}
	};

	for (m2 = m2_00_col; m2_end_cond; m2_diff) {
#ifdef DEBUG_BOUNDARY
		printf("S[%d][%d][sXX]=ZERO\n", m1_00_row, m2);
#endif

		for (i = 0; i < N_PAIR_HMM_STATES; i++) {
			S[m1_00_row][m2][i] = ZERO;
		}
	};

#ifdef DEBUG_BOUNDARY
	printf("all set up, running!\n\n");
#endif


	// the search itself
	if (dp->X_status == X_NOT_SET) {
#define X_STATUS X_NOT_SET
#include "DP_inner_wrapper_2.c"
#undef  X_STATUS

#if ALGORITHM != OPT_ACCUR
		dp->X_status = X_SET;
#endif
	} else if (dp->X_status == X_SET) {
#define X_STATUS X_SET
#include "DP_inner_wrapper_2.c"
#undef  X_STATUS
	} else if (dp->X_status == X_LIN2LOG) {
#define X_STATUS X_LIN2LOG
#include "DP_inner_wrapper_2.c"
#undef  X_STATUS

#if ALGORITHM != OPT_ACCUR
		dp->X_status = X_SET;
#endif
	};


#if (ALGORITHM != OPT_ACCUR) && (ALGO_SPACE == LINSPACE)
	// overflows are handled inside the inner loop
	if (dp->best_S < MIN_FLOAT) {
		dp->lin_overflow = 1;
	}
#endif

#if (ALGORITHM != OPT_ACCUR) && (ALGO_SPACE == LINSPACE)
	dp->best_S = log(dp->best_S);

#if ALGO_PATHS == SUM_PATHS
	dp->sum_S = log(dp->sum_S);
#endif
#endif

#ifdef DEBUG_SEARCH
	printf("best_S: %.2e @ (%d,%d)\n", dp->best_S, dp->best_m1, dp->best_m2);
#if ALGO_PATHS == SUM_PATHS
	printf("sum_S: %.2e\n", dp->sum_S);
#endif
#endif
}
