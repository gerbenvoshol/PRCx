/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

       DP_common.c: housekeeping functions and wrappers for DP_search.c

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
#include <string.h>
#include <ctype.h> /* Needed for tolower */
#include "prc.h"


// allocate the DPdata structure except S2, S_tr and the matches
//
DPdata* alloc_DPdata(PARAMS *p, HMM *hmm1, HMM *hmm2)
{
	DPdata *dp;

	calloc_1D_array(dp, DPdata, 1);

	dp->hmm1 = hmm1;
	dp->hmm2 = hmm2;
	dp->p    = p;

	malloc_3D_array(dp->S1, double,
	                hmm1->i->M + 4, hmm2->i->M + 4, N_PAIR_HMM_STATES);

	return dp;
}

// allocate X
//
void malloc_DPdata_X(DPdata *dp)
{
	if (dp->X == NULL)
		malloc_3D_array(dp->X, double,
		                dp->hmm1->i->M + 4, dp->hmm2->i->M + 4, N_X_STATES);
}

// the Waterman-Eggert algorithm: zero out those xMM states in the X matrix
// that lie on the trace, so that we can figure out what other matches there
// are
//
// return 1 if:
//
// - at least one of the values is already zero, meaning that we're re-crossing
//   an old path
//
// - the are *no* MM states anywhere on the trace; this isn't bad in itself,
//   but could lead to an infinite cycle, so we need to stop it
//
// otherwise, return 0
//
int zero_trace_X(DPdata *dp, MATCH *m)
{
	double test, zero = dp->lin_overflow ? LOG_ZERO : 0.0;
	int   i, n = 0;

	for (i = 1; i < m->n_trace - 1; i++) {
		if (m->trace[i][tr_state] != xMM) {
			continue;
		}

		n++;

		test = dp->X[ m->trace[i][tr_m1] ][ m->trace[i][tr_m2] ][xMM];
		dp->X[ m->trace[i][tr_m1] ][ m->trace[i][tr_m2] ][xMM] = zero;

		if (test == zero) {
			return 1;
		}
	};

	// n is the number of states that were zeroed
	if (n == 0) {
		return 1;
	}

	return 0;
}

// allocate S2
//
void malloc_DPdata_S2(DPdata *dp)
{
	if (dp->S2 == NULL)
		malloc_3D_array(dp->S2, double, dp->hmm1->i->M + 4, dp->hmm2->i->M + 4,
		                N_PAIR_HMM_STATES);
}

// allocate S_tr
//
void malloc_DPdata_S_tr(DPdata *dp)
{
	if (dp->S_tr == NULL)
		malloc_3D_array(dp->S_tr, int, dp->hmm1->i->M + 4, dp->hmm2->i->M + 4,
		                N_PAIR_HMM_STATES);
}

// add an empty structure at the start of the DPdata list of matches
//
// a pointer to the new structure is kept in DPdata.list
//
void add_match(DPdata *dp)
{
	MATCH *new;

	calloc_1D_array(new, MATCH, 1);

	new->hmm1i = dp->hmm1->i;
	new->hmm2i = dp->hmm2->i;

	new->next = dp->list;
	dp->list = new;
}

// fully free the DPdata structure, *except* the matches
//
// return a pointer to the start of the match list
//
MATCH* free_DPdata(DPdata *dp)
{
	MATCH *list = NULL;

	if (dp != NULL) {
		list = dp->list;

		if (dp->X_status == X_SET) {
			free_3D_array(dp->X);
		}

		free_3D_array(dp->S1);
		free_3D_array(dp->S2);
		free_3D_array(dp->S_tr);

		free(dp);
	};

	return list;
}

void print_DP(DPdata *dp, double ***S, int trace)
{
	// this will only work for SPACE9
	char *sXX[] = { "sMM", "sMI", "sIM", "sII",
	                "sMD", "sDM", "sDD", "sID", "sDI", "sBE"
	              };
	int m1, m2, i, tr;

	printf("\n");

	for (m1 = 0; m1 <= dp->hmm1->i->M + 3; m1++) {
		for (m2 = 0; m2 <= dp->hmm2->i->M + 3; m2++) {
			printf("m1=%d\tm2=%d\t", m1, m2);

			for (i = 0; i < N_PAIR_HMM_STATES; i++) {
				if (trace) {
					tr = dp->S_tr[m1][m2][i];
					if (tr < 0) {
						tr += 20;
					}
					printf("%s=%.1e(%s) ", sXX[i], S[m1][m2][i], sXX[tr]);
				} else {
					printf("%s=%.1e ", sXX[i], S[m1][m2][i]);
				};
			};

			printf("\n");
		};
	};

	printf("\n");
}

void write_Sxx_gnuplot(char *filename, DPdata *dp, double ***S, int state)
{
	FILE *file;
	int  m1, m2;

	open_file_or_die(file, filename, "w");

	for (m1 = 0; m1 <= dp->hmm1->i->M + 3; m1++) {
		for (m2 = 0; m2 <= dp->hmm2->i->M + 3; m2++)
			fprintf(file, "%d\t%d\t%.6f\n", m1, m2,
			        log(1e-20 + S[m1][m2][state]) / log(10));

		fprintf(file, "\n");
	};

	fclose(file);
}

// get the standard alignment string from the trace matrix
//
// i.e. something like:
//
// "MMMIMMM~MMM\n"
// "MMMMMMMDMMM"
//
// the correspondence with internal PRC states is as follows:
//
//          HMM1 HMM2
//          ~~~~ ~~~~
//   MM  ->  M    M
//   MI  ->  M    I
//   IM  ->  I    M
//   MD  ->  ~    D
//   ID  ->  ~    D 
//   DM  ->  D    ~
//   DI  ->  D    ~
//   DD  ->  D    D
//
void get_prc_alignment(MATCH *m)
{
	char *seq1, *seq2;
	long int n, i = 0;

	calloc_1D_array(seq1,         char, m->n_trace - 1);
	calloc_1D_array(seq2,         char, m->n_trace - 1);
	calloc_1D_array(m->alignment, char, m->n_trace * 2 - 1);

	for (n = 1; n < m->n_trace - 1; n++) {
		switch (m->trace[n][tr_state]) {
		case sMM:
			seq1[i] = 'M';
			seq2[i++] = 'M';
			break;
		case sMI:
			seq1[i] = 'M';
			seq2[i++] = 'I';
			break;
		case sIM:
			seq1[i] = 'I';
			seq2[i++] = 'M';
			break;
		case sMD:
			seq1[i] = '~';
			seq2[i++] = 'D';
			break;
		case sDM:
			seq1[i] = 'D';
			seq2[i++] = '~';
			break;
#if   PAIR_HMM_STATES == SPACE6
		case sDD:
			seq1[i] = 'D';
			seq2[i++] = 'D';
			break;
#elif PAIR_HMM_STATES == SPACE8
		case sDD:
			seq1[i] = 'D';
			seq2[i++] = 'D';
			break;
		case sID:
			seq1[i] = '~';
			seq2[i++] = 'D';
			break;
		case sDI:
			seq1[i] = 'D';
			seq2[i++] = '~';
			break;
#elif PAIR_HMM_STATES == SPACE9
		case sII:
			seq1[i] = 'I';
			seq2[i++] = 'I';
			break;
		case sDD:
			seq1[i] = 'D';
			seq2[i++] = 'D';
			break;
		case sID:
			seq1[i] = '~';
			seq2[i++] = 'D';
			break;
		case sDI:
			seq1[i] = 'D';
			seq2[i++] = '~';
			break;
#endif
		};
	};

	sprintf(m->alignment, "%s\n%s\n", seq1, seq2);

	free_unless_null(seq1);
	free_unless_null(seq2);
}

// get a SAM-like alignment string from the trace matrix
//
// i.e. something like:
//
// "MMmmMMDDMMMdddMMM----MMMMIIIMMMM"
//
// the correspondence with the standard alignment string above is as follows:
//
//   MM  ->  M
//   IM  ->  I
//   MI  ->  m
//   II  ->  i
//   ~D  ->  -
//   D~  ->  d
//   DD  ->  D
//   II  ->  i
//
// (for seq=1, i.e. when HMM1 is the "sequence")
//
void get_sam_alignment(MATCH *m, int seq)
{
	long int n, i = 0;
	char *str;

	calloc_1D_array(m->alignment, char, 1 + m->n_trace + m->n_trace / 60);
	str = m->alignment;

	if (seq == 1) {
		// HMM1 is the sequence
		for (n = 1; n < m->n_trace - 1; n++) {
			if (i % 61 == 60) {
				str[i++] = '\n';
			}

			switch (m->trace[n][tr_state]) {
			case sMM:
				str[i++] = 'M';
				break;
			case sMI:
				str[i++] = 'm';
				break;
			case sIM:
				str[i++] = 'I';
				break;
			case sMD:
				str[i++] = '-';
				break;
			case sDM:
				str[i++] = 'd';
				break;
#if   PAIR_HMM_STATES == SPACE6
			case sDD:
				str[i++] = 'D';
				break;
#elif PAIR_HMM_STATES == SPACE8
			case sDD:
				str[i++] = 'D';
				break;
			case sID:
				str[i++] = '-';
				break;
			case sDI:
				str[i++] = 'd';
				break;
#elif PAIR_HMM_STATES == SPACE9
			case sII:
				str[i++] = 'i';
				break;
			case sDD:
				str[i++] = 'D';
				break;
			case sID:
				str[i++] = '-';
				break;
			case sDI:
				str[i++] = 'd';
				break;
#endif
			};
		};
	} else if (seq == 2) {
		// HMM2 is the sequence => swap 'm' & 'I', 'd' & '-'
		for (n = 1; n < m->n_trace - 1; n++) {
			if (i % 61 == 60) {
				str[i++] = '\n';
			}

			switch (m->trace[n][tr_state]) {
			case sMM:
				str[i++] = 'M';
				break;
			case sMI:
				str[i++] = 'I';
				break;
			case sIM:
				str[i++] = 'm';
				break;
			case sMD:
				str[i++] = 'd';
				break;
			case sDM:
				str[i++] = '-';
				break;
#if   PAIR_HMM_STATES == SPACE6
			case sDD:
				str[i++] = 'D';
				break;
#elif PAIR_HMM_STATES == SPACE8
			case sDD:
				str[i++] = 'D';
				break;
			case sID:
				str[i++] = 'd';
				break;
			case sDI:
				str[i++] = '-';
				break;
#elif PAIR_HMM_STATES == SPACE9
			case sII:
				str[i++] = 'i';
				break;
			case sDD:
				str[i++] = 'D';
				break;
			case sID:
				str[i++] = 'd';
				break;
			case sDI:
				str[i++] = '-';
				break;
#endif
			};
		};
	} else {
		die2("The parameter 'seq' (=%d) must be either 1 or 2!", seq);
	};

	str[i] = '\n';

	//printf("\nAlignment:\n%s\n\n", str);
}

// set up hmm1,2->Smat,Sins based on the joint null
//
void setup_HMMs_for_dot2(HMM *hmm1, HMM *hmm2)
{
	double sqrt_null[20], sum = 0.0;
	int   i;

	for (i = 0; i < 20; i++) {
		sqrt_null[i] = hmm1->PinsJ[i] + hmm2->PinsJ[i];
		sum += sqrt_null[i];
	};

	for (i = 0; i < 20; i++) {
		sqrt_null[i] = sqrt(sqrt_null[i] / sum);
	}

	scorify_HMM_M_I(hmm1, sqrt_null);
	scorify_HMM_M_I(hmm2, sqrt_null);
}

// set alignment modes for the two HMMs
//
// the B->Mx and Mx->E modes are set to LOCAL or GLOBAL independently
//
void set_HMM_align_modes(DPdata* dp, REGION *reg,
                         int BM1, int ME1, int BM2, int ME2)
{
	set_HMM_align_mode(dp->hmm1, tBM, BM1, reg->start1, reg->end1);
	set_HMM_align_mode(dp->hmm1, tME, ME1, reg->start1, reg->end1);
	set_HMM_align_mode(dp->hmm2, tBM, BM2, reg->start2, reg->end2);
	set_HMM_align_mode(dp->hmm2, tME, ME2, reg->start2, reg->end2);
}

// a wrapper for doing a simple run
//
// returns:
//
// * 1 if all goes well
//
// * 666 if there is a LINSPACE overflow
//
int DP_wrapper(DPdata *dp)
{
	set_HMM_align_modes
	(dp, &dp->region,
	 dp->p->align_mode1, dp->p->align_mode1,
	 dp->p->align_mode2, dp->p->align_mode2);

	if (dp->p->algorithm == PARAM_FORW_BACK) {
		if (dp->lin_overflow) {
			search_LOGSPACE_FORWARD(dp, &dp->region);
		} else {
			//write_Sxx_gnuplot("/ss0/mm238/back.dat",  dp, dp->S2, sMM);
			//print_DP(dp, dp->S1, 0);

			search_LINSPACE_FORWARD(dp, &dp->region);

			if (dp->lin_overflow) {
				return 666;
			}
		};
	} else if (dp->p->algorithm == PARAM_VITERBI) {
		if (dp->lin_overflow) {
			search_LOGSPACE_VITERBI(dp, &dp->region);
			//print_DP(dp, dp->S1, 1);
		} else {
			search_LINSPACE_VITERBI(dp, &dp->region);
			if (dp->lin_overflow) {
				return 666;
			}
		};
	} else {
		die2("Unknown algorithm #%d!", dp->p->algorithm);
	};

	return 1;
}

// finish a forward/backward run started by DP_wrapper() above
//
// returns:
//
// * 0 if DP_wrapper() ran Viterbi (in which case this does nothing)
//
// * 1 if DP_wrapper() ran forward and there are no problems finishing off
//
// * 666 if DP_wrapper() ran forward in LINSPACE without overflow and one of
//   the runs in this routine overflows; also sets dp->lin_overflow
//
// NEED TO EXPLAIN WHAT I'M DOING HERE!
//
//#define DEBUG_FINISH
int DP_finish_forw_back(DPdata* dp)
{
	REGION match, search_region;
	int    i;

	// N.B. dp->algorithm is set to the last algorithm run
	if (dp->algorithm != FORWARD) {
		return 0;
	}

	copy_region(&dp->region, &match);

#define LOTS 5
	//
	// This is a completely ad hoc number; ideally the for loop should be an
	// infinite loop, but ... infinite loops are nasty, plus I've seen cases
	// where the behaviour is oscillatory and there's no convergence
	//
	// [though this may have been due to a bug in BACKWARD, fixed in 1.5.0]
	//
	// 5 should be more than enough for convergence
	//
	for (i = 0; i < LOTS; i++) {
		// BACKWARD loop

		match.end1 = dp->best_m1;
		match.end2 = dp->best_m2;

#ifdef DEBUG_FINISH
		printf("iteration %d end: %d,%d\n", i, match.end1, match.end2);
#endif

		// the start is as before, but the end is at (best_m1, best_m2)
		//
		set_region(&search_region,
		           dp->region.start1, match.end1,
		           dp->region.start2, match.end2);

		set_HMM_align_modes
		(dp, &search_region,
		 dp->p->align_mode1, dp->p->align_mode1,
		 dp->p->align_mode2, dp->p->align_mode2);

		if (dp->lin_overflow) {
			search_LOGSPACE_BACKWARD(dp, &search_region);

			//printf("backward matrix\n");
			//print_DP(dp, dp->S2, 0);
		} else {
			search_LINSPACE_BACKWARD(dp, &search_region);

			if (dp->lin_overflow) {
				return 666;
			}
		};

#ifdef DEBUG_FINISH
		printf("iteration %d BACKWARD scores: sum=%f best=%f\n",
		       i, dp->sum_S, dp->best_S);
#endif

		//printf("backward matrix\n");
		//print_DP(dp, dp->S2, 0);
		//write_Sxx_gnuplot("sMM.dat",  dp, dp->S1, sMM);
		//write_Sxx_gnuplot("backward.dat", dp, dp->S2, sMM);

		// convergence
		if ((dp->best_m1 == match.start1) && (dp->best_m2 == match.start2)) {
			break;
		}

		match.start1 = dp->best_m1;
		match.start2 = dp->best_m2;

#ifdef DEBUG_FINISH
		printf("iteration %d start: %d,%d\n", i, match.start1, match.start2);
#endif

		set_region(&search_region,
		           match.start1, dp->region.end1,
		           match.start2, dp->region.end2);

		set_HMM_align_modes
		(dp, &search_region,
		 dp->p->align_mode1, dp->p->align_mode1,
		 dp->p->align_mode2, dp->p->align_mode2);

		if (dp->lin_overflow) {
			search_LOGSPACE_FORWARD(dp, &search_region);
		} else {
			search_LINSPACE_FORWARD(dp, &search_region);

			if (dp->lin_overflow) {
				return 666;
			}
		};

#ifdef DEBUG_FINISH
		printf("iteration %d FORWARD scores: sum=%f best=%f\n",
		       i, dp->sum_S, dp->best_S);
#endif

		// convergence
		if ((dp->best_m1 == match.end1) && (dp->best_m2 == match.end2)) {
			break;
		}
	};

	//if( i==LOTS )
	//  { die1("Uh-oh, i=LOTS. This really shouldn't happen ..."); };

	/*
	printf("opt_accur matrix\n");
	if( dp->algorithm == FORWARD )
	  print_DP(dp, dp->S2, 0);
	else if( dp->algorithm == BACKWARD )
	  print_DP(dp, dp->S1, 0);
	*/

#ifdef DEBUG_FINISH
	printf("Done!\n");
	printf("Match: %d-%d, %d-%d\n\n",
	       match.start1, match.end1,
	       match.start2, match.end2);
#endif

	copy_region(&match, &dp->list->in_out);

	return 1;
}

// do a traceback for VITERBI or OPT_ACCUR
//
// returns:
//
// * 1 if all goes well
//
// * 666 if BACKWARD or OPT_ACCUR overflow in LINSPACE
//
int traceback_wrapper(DPdata *dp)
{
	// N.B. dp->algorithm is set to the last algorithm run
	if (dp->algorithm == VITERBI) {
		traceback_VITERBI(dp, dp->best_m1, dp->best_m2);
	} else if (dp->algorithm == FORWARD) {
		if (DP_finish_forw_back(dp) == 666) {
			return 666;
		}

		if (dp->lin_overflow) {
			search_LOGSPACE_OPT_ACCUR(dp, &dp->list->in_out);
		} else {
			search_LINSPACE_OPT_ACCUR(dp, &dp->list->in_out);

			if (dp->lin_overflow) {
				return 666;
			}
		};

		traceback_OPT_ACCUR(dp, dp->list->in_out.end1, dp->list->in_out.end2);
	} else {
		die2("Unknown algorithm #%d!", dp->algorithm);
	};

	return 1;
}

//#define DEBUG_RUN
void log_run(DPdata *dp)
{
#ifdef DEBUG_RUN
	printf("Entering log_run() ...\n");
#endif

	dp->X_status = X_LIN2LOG;
	DP_wrapper(dp);

	if (dp->p->algorithm == PARAM_FORW_BACK) {
		dp->list->sum = dp->sum_S;
	}

	traceback_wrapper(dp);
};

void lin_run(DPdata *dp)
{
	REGION *r;
	int    length;
	double  log_length;

#ifdef DEBUG_RUN
	printf("Entering lin_run() ...\n");
#endif

	add_match(dp);


	if (DP_wrapper(dp) == 666) {
		log_run(dp);
	} else {
		if (dp->p->algorithm == PARAM_FORW_BACK) {
			dp->list->sum = dp->sum_S;
		}
		/* If an overflow/underflow was detected do a logarithmic run */
		if (traceback_wrapper(dp) == 666) {
			log_run(dp);
		}
	};

#ifdef DEBUG_RUN
	printf("Did the first DP_wrapper ...\n");
#endif

	if (dp->p->algorithm == PARAM_VITERBI) {
		dp->list->simple = dp->best_S;
	} else if (dp->p->algorithm == PARAM_FORW_BACK) {
		r = &(dp->list->proper);
		length = r->end1 - r->start1 + 1 + r->end2 - r->start2 + 1;
		log_length = (double) log((0.5 * ((double)length)));

		if (dp->lin_overflow) {
			dp->list->simple = dp->best_S - log_length;
		} else {
			dp->list->simple = log(dp->best_S) - log_length;
		}


#ifdef DEBUG_RUN
		printf("lin_run: length=%d best_S=%e (overfl=%d) simple=%f\n",
		       length, dp->best_S, dp->lin_overflow, dp->list->simple);
#endif

	};
};

/* TODO: Place in correct file HMM.c? */
char *get_ems_consensus(HMM *hmm)
{
	int i, n;
	/* The amino acid alphabet */
	const char AA[20] = "ACDEFGHIKLMNPQRSTVWY";
	/* The default background frequencies */
	/*                         A         C         D         E        F         G         H         I         K        L         M         N         P         Q         R         S         T         V         W         Y      */
	//const double nule[20] = {0.078785, 0.015156, 0.053514, 0.06685, 0.039694, 0.069496, 0.022925, 0.059009, 0.05946, 0.096393, 0.023766, 0.041437, 0.048297, 0.039557, 0.054111, 0.068349, 0.054074, 0.067362, 0.011415, 0.030418};
	char *consensus = calloc(hmm->i->M + 1, sizeof(char));
	double max_value;
	int max_index;

	/* Go through the emission scores and find the AA with the highest score */
	for (i = 2; i < hmm->i->M + 2; i++) {
		/* Find the AA with the highest emsission */
		max_value = 0.0;
		max_index = 0;
		for (n = 0; n < 20; n++) {
			/* Store the value of the highest emission and the corresponding index */
			if (hmm->Pmat[i][n] > max_value) {
				max_value = hmm->Pmat[i][n];
				max_index = n;
			}
		}

		/* Construct the consensus line, return an upper case letter if the probability is 
		 * higher then expected from the background distribution and a lower case letter otherwise
		 */
		if (max_value > 0.5) {
			consensus[i - 2] = AA[max_index];
		} else {
			consensus[i - 2] = tolower(AA[max_index]);
		}
	}

	/* Done */
	return consensus;
}

// get the standard alignment string from the trace matrix
//
// i.e. something like:
//
// "MMMIMMM~MMM\n"
// "MMMMMMMDMMM"
//
void get_prc2_alignment(MATCH *m, HMM *hmm1, HMM *hmm2)
{
	/* Get the consensus string of the HMM */
	char *ems_cons_hmm1 = NULL;
	char *ems_cons_hmm2 = NULL;
	ems_cons_hmm1 = get_ems_consensus(hmm1);
	ems_cons_hmm2 = get_ems_consensus(hmm2);

	char *seq1, *seq2;
	long int n, i = 0, l1 = 0, l2 = 0;

	calloc_1D_array(seq1,         char, m->n_trace - 1);
	calloc_1D_array(seq2,         char, m->n_trace - 1);
	calloc_1D_array(m->alignment, char, m->n_trace * 2 - 1);

	for (n = 1; n < m->n_trace - 1; n++) {
		switch (m->trace[n][tr_state]) {
		case sMM:
			seq1[i] = ems_cons_hmm1[l1 + m->proper.start1 - 2];
			seq2[i++] = ems_cons_hmm2[l2 + m->proper.start2 - 2];
			l1++;
			l2++;
			break;
		case sMI:
			seq1[i] = ems_cons_hmm1[l1 + m->proper.start1 - 2];
			seq2[i++] = '.';
			l1++;
			break;
		case sIM:
			seq1[i] = '.';
			seq2[i++] = ems_cons_hmm2[l2 + m->proper.start2 - 2];
			l2++;
			break;
		case sMD:
			seq1[i] = '-';
			seq2[i++] = ems_cons_hmm2[l2 + m->proper.start2 - 2];
			l2++;
			break;
		case sDM:
			seq1[i] = ems_cons_hmm1[l1 + m->proper.start1 - 2];
			seq2[i++] = '-';
			l1++;
			break;
#if   PAIR_HMM_STATES == SPACE6
		case sDD:
			seq1[i] = '-';
			seq2[i++] = '-';
			break;
#elif PAIR_HMM_STATES == SPACE8
		case sDD:
			seq1[i] = '-';
			seq2[i++] = '-';
			break;
		case sID:
			seq1[i] = '.';
			seq2[i++] = '-';
			break;
		case sDI:
			seq1[i] = '-';
			seq2[i++] = '.';
			break;
#elif PAIR_HMM_STATES == SPACE9
		case sII:
			seq1[i] = '.';
			seq2[i++] = '.';
			break;
		case sDD:
			seq1[i] = '-';
			seq2[i++] = '-';
			l1++;
			l2++;
			break;
		case sID:
			seq1[i] = '.';
			seq2[i++] = '-';
			l2++;
			break;
		case sDI:
			seq1[i] = '-';
			seq2[i++] = '.';
			l1++;
			break;
#endif
		};
	};

	sprintf(m->alignment, "%s\n%s\n", seq1, seq2);

	/* Cleanup */
	free(ems_cons_hmm1);
	free(ems_cons_hmm2);

	free_unless_null(seq1);
	free_unless_null(seq2);
}


// do a full run of two HMMs against each other, including the reverse null
//
// the single_match parameter determines whether to get just the best match
// (for single_match=1), e.g. in order to estimate E-value statistics for all
// models in a library, or whether to get all hits satisfying the various
// thresholds
//
// if the thresholds are used, the routine may, and probably will, return some
// hits that don't satisfy them
//
// the do_alignments parameter tells the routine whether or not it should call
// the alignment routines
//
// returns a pointer to the DPdata list of matches
//
MATCH* run_two_HMMs(PARAMS *p, HMM *hmm1, HMM *rev1, HMM *hmm2,
                    int single_match, int do_alignments)
{
	DPdata *dp = NULL, *dp_rev = NULL;
	MATCH  *m;
	long   int hit_i = 0;
	double  reverse = 0.0;

#ifdef DEBUG_RUN
	printf("Entering run_two_HMMs() ...\n");
#endif

	// do the reverse null
#ifdef DEBUG_RUN
	printf("Doing the reverse null ...\n");
#endif
	if (p->MM_function == DOT2) {
		setup_HMMs_for_dot2(rev1, hmm2);
	}

	dp_rev = alloc_DPdata(p, rev1, hmm2);
	set_region(&dp_rev->region, 1, hmm1->i->M + 2, 1, hmm2->i->M + 2);

	lin_run(dp_rev);
	reverse = dp_rev->list->simple;

	free_matches(dp_rev->list);
	free_DPdata(dp_rev);

#ifdef DEBUG_RUN
	printf("Done the reverse ...\n");
#endif


	if (p->MM_function == DOT2) {
		setup_HMMs_for_dot2(hmm1, hmm2);
	}

	dp = alloc_DPdata(p, hmm1, hmm2);
	set_region(&dp->region, 1, hmm1->i->M + 2, 1, hmm2->i->M + 2);

	do {
#ifdef DEBUG_RUN
		printf("Starting do loop ...\n");
#endif

		// the forward match
		lin_run(dp);

		m = dp->list;
		m->match_n = ++hit_i;
		m->reverse = m->simple - reverse;

#ifdef DEBUG_RUN
		printf("Done traceback_wrapper() ...\n");
		printf("in_out: %d-%d, %d-%d\n",
		       m->in_out.start1, m->in_out.end1,
		       m->in_out.start2, m->in_out.end2);
		printf("proper: %d-%d, %d-%d\n",
		       m->proper.start1, m->proper.end1,
		       m->proper.start2, m->proper.end2);
#endif

		if (do_alignments) {
			switch (p->align_style) {
			case ALIGN_PRC1:
				get_prc_alignment(m);
				break;
			case ALIGN_PRC2:
				get_prc2_alignment(m, hmm1, hmm2);
				break;
			case ALIGN_SAM1:
				get_sam_alignment(m, 1);
				break;
			case ALIGN_SAM2:
				get_sam_alignment(m, 2);
				break;
			case ALIGN_NONE:
				break;
			default:
				die2("Uknown alignment style #%d", p->align_style);
			};
		};

#ifdef DEBUG_RUN
		printf("Done get_XXX_alignment() ...\n");
#endif

		if (single_match) {
			break;
		}

		if (p->E_values) {
			m->E_value = E_value(p->n_unrel, m->reverse, p->p);
		}

		// if zero_trace_X returns 1, it means that this path crossed an
		// earlier one (or there are no MM states), so this sounds like a
		// good time to stop
		if (zero_trace_X(dp, m)) {
			// if hit_i==1, then there are no MM states but we need to return
			// *something*, so don't delete the match
			if (hit_i > 1) {
				dp->list = m->next;
				free_match(m);
			};

			break;
		};

#ifdef DEBUG_RUN
		printf("Match model=%s stop=%f:\n"
		       "simple=%f reverse=%f E=%e, 1:%d-%d, 2:%d-%d\n",
		       m->hmm2i->name,
		       p->stop,
		       m->simple, m->reverse, m->E_value,
		       m->proper.start1, m->proper.end1,
		       m->proper.start2, m->proper.end2);

		if (p->align_style != ALIGN_NONE) {
			printf("%s\n", m->alignment);
		}
#endif
	} while (WEAK_MATCH(p, m) && (hit_i < p->max_hits));


#ifdef DEBUG_RUN
	printf("Leaving run_two_HMMs() ...\n\n");
#endif

	// free_DPdata returns the list of matches
	return free_DPdata(dp);
}

// do a simple run of two HMMs against each other, without the reverse null
//
// the single_match parameter determines whether to get just the best match
// (for single_match=1), e.g. in order to estimate E-value statistics for all
// models in a library, or whether to get all hits satisfying the various
// thresholds
//
// if the thresholds are used, the routine may, and probably will, return some
// hits that don't satisfy them
//
// the do_alignments parameter tells the routine whether or not it should call
// the alignment routines
//
// returns a pointer to the DPdata list of matches
//
MATCH* simple_run_two_HMMs(PARAMS *p, HMM *hmm1, HMM *hmm2,
                    int single_match, int do_alignments)
{
	DPdata *dp = NULL;
	MATCH  *m;
	long   int hit_i = 0;

#ifdef DEBUG_RUN
	printf("Entering simple_run_two_HMMs() ...\n");
#endif

	if (p->MM_function == DOT2) {
		setup_HMMs_for_dot2(hmm1, hmm2);
	}

	dp = alloc_DPdata(p, hmm1, hmm2);
	set_region(&dp->region, 1, hmm1->i->M + 2, 1, hmm2->i->M + 2);

	do {
#ifdef DEBUG_RUN
		printf("Starting do loop ...\n");
#endif

		// the forward match
		lin_run(dp);

		m = dp->list;
		m->match_n = ++hit_i;
		// We did not do a run against the reverse HMM!
		m->reverse = 0.0;

#ifdef DEBUG_RUN
		printf("Done traceback_wrapper() ...\n");
		printf("in_out: %d-%d, %d-%d\n",
		       m->in_out.start1, m->in_out.end1,
		       m->in_out.start2, m->in_out.end2);
		printf("proper: %d-%d, %d-%d\n",
		       m->proper.start1, m->proper.end1,
		       m->proper.start2, m->proper.end2);
#endif

		if (do_alignments) {
			switch (p->align_style) {
			case ALIGN_PRC1:
				get_prc_alignment(m);
				break;
			case ALIGN_PRC2:
				get_prc2_alignment(m, hmm1, hmm2);
				break;
			case ALIGN_SAM1:
				get_sam_alignment(m, 1);
				break;
			case ALIGN_SAM2:
				get_sam_alignment(m, 2);
				break;
			case ALIGN_NONE:
				break;
			default:
				die2("Uknown alignment style #%d", p->align_style);
			};
		};

#ifdef DEBUG_RUN
		printf("Done get_XXX_alignment() ...\n");
#endif

		if (single_match) {
			break;
		}

		if (p->E_values) {
			// Since we did not run the reverse, we cannot calculate a meaningful E-value
			m->E_value = 0.0;
		}

		// if zero_trace_X returns 1, it means that this path crossed an
		// earlier one (or there are no MM states), so this sounds like a
		// good time to stop
		if (zero_trace_X(dp, m)) {
			// if hit_i==1, then there are no MM states but we need to return
			// *something*, so don't delete the match
			if (hit_i > 1) {
				dp->list = m->next;
				free_match(m);
			};

			break;
		};

#ifdef DEBUG_RUN
		printf("Match model=%s stop=%f:\n"
		       "simple=%f reverse=%f E=%e, 1:%d-%d, 2:%d-%d\n",
		       m->hmm2i->name,
		       p->stop,
		       m->simple, m->reverse, m->E_value,
		       m->proper.start1, m->proper.end1,
		       m->proper.start2, m->proper.end2);

		if (p->align_style != ALIGN_NONE) {
			printf("%s\n", m->alignment);
		}
#endif
	} while (WEAK_MATCH(p, m) && (hit_i < p->max_hits));


#ifdef DEBUG_RUN
	printf("Leaving simple_run_two_HMMs() ...\n\n");
#endif

	// free_DPdata returns the list of matches
	return free_DPdata(dp);
}

