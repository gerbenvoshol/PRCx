/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

		   DP_traceback.c: traceback for DP_search.c

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prc.h"


#if   ALGORITHM == VITERBI
#define traceback  traceback_VITERBI
#elif ALGORITHM == OPT_ACCUR
#define traceback  traceback_OPT_ACCUR
#endif

//#define DEBUG_TRACEBACK
void traceback(DPdata *dp, int start1, int start2)
{
	MATCH *match = dp->list;
	int   n, m, m1, m2, state, prev_state, **trace;

	// this is only going to work for SPACE9, edit as necessary
#ifdef DEBUG_TRACEBACK
	char *sXX[]   = { "sMM", "sMI", "sIM", "sII",
	                  "sMD", "sDM", "sDD", "sID", "sDI", "sBE"
	                };
#endif

	m1      = start1;
	m2      = start2;
	n       = 0;
	state   = sMM;

	malloc_DPdata_S_tr(dp);
	malloc_2D_array(trace, int, match->hmm1i->M + match->hmm2i->M + 2, 3);

#ifdef DEBUG_TRACEBACK
	printf("\n\nTRACE:\n\n");
	printf("  m1 = %d, m2 = %d\n\n", m1, m2);
#endif

#if   ALGORITHM == VITERBI
#define cond (state != sBE)
#elif ALGORITHM == OPT_ACCUR
#define cond ((m1>=match->in_out.start1) && (m2>=match->in_out.start2))
#endif
	while (cond) {
		trace[n][tr_m1]    = m1;
		trace[n][tr_m2]    = m2;
		trace[n][tr_state] = state;

		prev_state = dp->S_tr[m1][m2][state];

#if ALGORITHM == VITERBI
		if (prev_state < 0) {
			if (n == 0) {
				prev_state += 20;
			} else {
				prev_state = sBE;
			}
		};
#endif

#ifdef DEBUG_TRACEBACK
		printf("%d:\t%s -> %s(%d,%d)\t = %.6f\n",
		       n,
		       sXX[prev_state],
		       sXX[state], m1, m2,
		       dp->S1[m1][m2][state]);
#endif

		switch (state) {
		case sMM:
			m1--;
			m2--;
			break;
		case sMD:
			m2--;
			break;
		case sDM:
			m1--;
			break;
		case sMI:
			m1--;
			break;
		case sIM:
			m2--;
			break;
#if   PAIR_HMM_STATES == SPACE6
		case sDD:
			m1--;
			m2--;
			break;
#elif PAIR_HMM_STATES == SPACE8
		case sDD:
			m1--;
			m2--;
			break;
		case sID:
			m2--;
			break;
		case sDI:
			m1--;
			break;
#elif PAIR_HMM_STATES == SPACE9
		case sII:
			break;
		case sDD:
			m1--;
			m2--;
			break;
		case sID:
			m2--;
			break;
		case sDI:
			m1--;
			break;
#endif
		};

		state = prev_state;
		n++;
	};

#ifdef DEBUG_TRACEBACK
	printf("trace done ...\n");
#endif

	// save the reverse of **trace
	match->n_trace = n;
	malloc_2D_array(match->trace, int, n, 3);

	for (m = 0; m < n; m++) {
		match->trace[n - 1 - m][tr_m1]    = trace[m][tr_m1];
		match->trace[n - 1 - m][tr_m2]    = trace[m][tr_m2];
		match->trace[n - 1 - m][tr_state] = trace[m][tr_state];
	};

	free_2D_array(trace);

	set_region(&match->proper,
	           match->trace[1][tr_m1], match->trace[match->n_trace - 2][tr_m1],
	           match->trace[1][tr_m2], match->trace[match->n_trace - 2][tr_m2]);

	set_region(&match->in_out,
	           match->trace[0][tr_m1], match->trace[match->n_trace - 1][tr_m1],
	           match->trace[0][tr_m2], match->trace[match->n_trace - 1][tr_m2]);

#ifdef DEBUG_TRACEBACK
	printf("all done ...\n");
#endif
}
