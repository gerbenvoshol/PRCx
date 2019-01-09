/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   PRCx, the profile comparer eXtended, version 1.0.0

		   stats.c: E-value fitting and calculation

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
#include <string.h>
#include <math.h>
#include "prc.h"

//#define DEBUG
//#define STAND_ALONE

#ifdef STAND_ALONE
//
// parse a PRC output file
//
void read_file(char *name,  long int *n, double **xs)
{
	FILE   *file;
	int    aux_int, hit_no;
	double aux_double, reverse;
	long   int aux_long, i = 0;
	char   line[500], buffer[500];

	if (!(file = fopen(name, "r"))) {
		char error[500];
		sprintf(error, "Error opening file '%s'", name);
		perror(error);
		exit(1);
	};

	while (fgets(line, 500, file)) {
		if (sscanf(line, "# Lib size n : %ld", &aux_long)) {
			*n = aux_long;

#ifdef DEBUG
			printf("Got n=%ld\n\n", *n);
#endif
			malloc_1D_array(*xs,    double, *n);

			continue;
		};

		if (line[0] == '#') {
			continue;
		}

		sscanf(line,

		       // 1=hmm1 2=start1 3=end1 4=length1 5=hit_no
		       //1   2   3   4   5
		       "%s\t%d\t%d\t%d\t%d\t"

		       // 6=hmm2 7=start2 8=end2 9=length2
		       //6   7   8   9
		       "%s\t%d\t%d\t%d\t"

		       // 10=co-emis 11=simple  12=reverse  13=E-value
		       //10   11   12  13
		       "%lf\t%lf\t%lf\t%lf\n",

		       buffer, &aux_int, &aux_int, &aux_int, &hit_no,
		       buffer, &aux_int, &aux_int, &aux_int,
		       &aux_double, &aux_double, &reverse, &aux_double);

#ifdef DEBUG
		//printf("%ld:\t%d\t%f\n", i, hit_no, reverse);
#endif

		if (hit_no > 1) {
			continue;
		}

		(*xs)[i]    = reverse;
		i++;
	};

	*n = i;
#ifdef DEBUG
	printf("\nSet n=%d\n\n", *n);
#endif

	fclose(file);
}
#endif

// Return the gradient of the log-likelihood function
//
// The derivation is given in my thesis in the PRC chapter, the section on
// Statistical significance. The final equations are:
//
//
//    pL/pl = n_bot * x_bot / [ 1 + exp(l*x_bot+k) ]  +  n_fit / l
//
//            + Sum(i=1..n_fit) x_i * tanh[ -0.5*(l*x_i-k) ]
//
//    pL/pk = n_bot / [ 1 + exp(l*x_bot+k) ]  +
//
//            + Sum(i=1..n_fit) tanh[ -0.5*(l*x_i-k) ]
//
//
// where the variables and parameters are:
//
//      L  ...  log-likelihood function
//      l  ...  lambda, p[0]
//      k  ...  kappa,  p[1]
//  n_top  ...  significant matches, ignore
//  n_fit  ...  matches to fit
//  n_bot  ...  ignore n_bot bottom matches (possibly poor statistics)
//  x_bot  ...  score of the bottom fitted match, xs[n_top+n_fit-1]
//  pL/pl  ...  gr[0], partial dL / partial dlambda
//  pL/pk  ...  gr[1], partial dL / partial dkappa
//
//
void grad(long int n_top, long int n_fit, long int n_bot,
          double *xs, double *p, double *gr)
{
	double xi, w, th, fr;
	long   int i;

	gr[0] = gr[1] = 0.0;

	//printf("n_top=%ld, n_fit=%ld, n_bot=%ld\n", n_top, n_fit, n_bot);

	for (i = n_top; i < n_top + n_fit; i++) {
		xi    = xs[i];
		w     = p[0] * xi + p[1];
		th    = tanh(-0.5 * w);

		//printf("xi=%f, w=%f, th=%f\n", xi, w, th);

		gr[0] += th * xi;
		gr[1] += th;
	};

	xi = xs[n_top + n_fit - 1];
	w  = p[0] * xi + p[1];
	fr = ((double)n_bot) / (1 + exp(w));

	gr[0] += fr * xi + ((double)n_fit) / p[0];
	gr[1] += fr;

#ifdef DEBUG
	printf("the gradient at (%f,%f) is (%f,%f)\n",
	       p[0], p[1], gr[0], gr[1]);
#endif
}

// find the zero of _d.grad(f)_ along _d_ from _p_ using the secand method
//
void zero_dot(long int n_top, long int n_fit, long int n_bot,
              double *xs, double *p,  double *d,  double *gr)
{
	double alpha = 1.0e-7, max_alpha, dot, dot_old;
	long   int i;

#ifdef DEBUG
	printf("searching along (%f,%f) from (%f,%f) ... \n",
	       d[0], d[1], p[0], p[1]);
#endif

	// make sure that |alpha * d|^2 < 0.1
	max_alpha = 0.1 / (d[0] * d[0] + d[1] * d[1]);
	while (alpha * alpha > max_alpha) {
		alpha /= 10.0;
	}

	p[0] -= alpha * d[0];
	p[1] -= alpha * d[1];

	grad(n_top, n_fit, n_bot, xs, p, gr);
	dot_old = d[0] * gr[0] + d[1] * gr[1];
#ifdef DEBUG
	printf("dot(%f,%f) = %f\n", p[0], p[1], dot_old);
#endif

	for (i = 0; i < 100; i++) {
		while (alpha * alpha > 0.1 * max_alpha) {
			alpha /= 10.0;
		}

		p[0] += alpha * d[0];
		p[1] += alpha * d[1];
		grad(n_top, n_fit, n_bot, xs, p, gr);
		dot = d[0] * gr[0] + d[1] * gr[1];
#ifdef DEBUG
		printf("dot(%f,%f) = %f\n", p[0], p[1], dot);
#endif

		if (fabs(dot) < 1.0e-10) {
#ifdef DEBUG
			printf("dot < 1e-10, search over\n\n");
#endif
			break;
		};

		if (dot_old == dot) {
#ifdef DEBUG
			printf("dot won't budge, search over\n\n");
#endif
			break;
		};

		alpha *= dot / (dot_old - dot);
		dot_old = dot;
	};

#ifdef DEBUG
	if (i == 100) {
		printf("i=100, the method failed to converge\n\n");
	}
#endif
}

// conjugate gradient minimization of the log-likelihood function
//
void conj_grad(long int n_top, long int n_fit, long int n_bot,
               double *xs, double *p)
{
	double dot_old, dot_mix, dot_new, alpha = 0.0;
	double old_gr[2], gr[2], d[2] = {0.0, 0.0};
	long   int i;

	grad(n_top, n_fit, n_bot, xs, p, old_gr);
	dot_old = old_gr[0] * old_gr[0] + old_gr[1] * old_gr[1];

	for (i = 1; i <= 50; i++) {
		d[0] = alpha * d[0] - old_gr[0];
		d[1] = alpha * d[1] - old_gr[1];

		zero_dot(n_top, n_fit, n_bot, xs, p, d, gr);

		dot_new = gr[0] * gr[0] + gr[1] * gr[1];
		dot_mix = gr[0] * old_gr[0] + gr[1] * old_gr[1];
		alpha = (dot_new - dot_mix) / dot_old;

#ifdef DEBUG
		printf("at (%f,%f)\n", p[0], p[1]);
		printf(" ..  grad f = (%f,%f)\n", gr[0], gr[1]);
		printf(" .. dot_new = %e\n", dot_new);
		printf("\n");
#endif

		if (fabs(dot_new) < 1.0e-6) {
			break;
		}

		old_gr[0] = gr[0];
		old_gr[1] = gr[1];
		dot_old = dot_new;

		if (alpha < 0.0) {
#ifdef DEBUG
			printf("alpha less than zero! restarting ...\n\n");
#endif
			alpha = 0.0;
		} else if (i % 5 == 0) {
#ifdef DEBUG
			printf("five iterations, restarting ... \n\n");
#endif
			alpha = 0.0;
		};
	};
}

// calculate the E-value for:
//
//   n_unrel ...  number of unrelated targets
//   S       ...  the query-target score
//   p       ...  {lambda, kappa} = parameters of the distribution
//
//
// the Karplus function for E-values based on the reverse null is:
//
//                           1
//     P(S>x) = ---------------------------
//               1 + exp(lambda*x + kappa)
//
// and
//
//       E(S) = n_unrel * P(S>x)
//
double E_value(long int n_unrel, double S, double *p)
{
	return ((double)n_unrel) / (1.0 + exp(p[0] * S + p[1]));
}


// for an E-value argument, return a nicely formatted string with a sane number
// of significant figures
//
char* nice_string(double E)
{
	char *string;

	calloc_1D_array(string, char, 10);

	if (E < 1.0e-4) {
		sprintf(string, "%.1e", E);
	} else if (E < 0.000995) {
		sprintf(string, "%.5f", E);
	} else if (E < 0.00995) {
		sprintf(string, "%.4f", E);
	} else if (E < 0.0995) {
		sprintf(string, "%.3f", E);
	} else if (E < 0.995) {
		sprintf(string, "%.2f", E);
	} else if (E < 9.95) {
		sprintf(string, "%.1f", E);
	} else if (E < 999.5) {
		sprintf(string, "%.0f", E);
	} else {
		sprintf(string, "%.1e", E);
	}

	return string;
}

// remove hits with E-value < 1.0, and bottom 10% of hits with E-value >= 1.0,
// to avoid biasing the statistics by non-random true matches at the top and
// poorly fitting weak matches at the bottom
//
//    n_top = number of hits with E-value < 1.0
//    n_fit = number of hits actually fitted
//    n_bot = number of bottom 10% unrelated hits
//
// together,
//
//    n_top + n_fit + n_bot = n
//
// and
//
//    n_bot / (n_fit + n_bot) >= 10%
//
// for E_value() above,
//
//    n_unrel = 1 + n_fit + n_bot
//
//
// NB it is assumed that *xs is sorted
//
#define FR_FIT 0.9
void censor(long int n, double *xs,
            long int *n_top, long int *n_fit, long int *n_bot,
            double *p)
{
	double E;
	long   int i, new_n_top = 0;

	for (i = 0; i < n; i++) {
		E = E_value(1 + *n_fit + *n_bot, xs[i], p);

#ifdef DEBUG
		if (E < 1.0) {
			printf("i=%ld: Dropping E=%e for x=%f\n", i, E, xs[i]);
		}
#endif

		if (E < 1.0) {
			new_n_top++;
		} else {
			break;
		}
	};

	*n_top = new_n_top;
	*n_fit = (long int)(FR_FIT * ((double)(n - *n_top)));
	*n_bot = n - *n_fit - *n_top;

#ifdef DEBUG
	printf("\n");
	printf("... n_top = %ld\n", *n_top);
	printf("... n_fit = %ld\n", *n_fit);
	printf("... n_bot = %ld\n", *n_bot);
	printf("\n");
#endif
}

// called by iterate_censor_fit() below to qsort *xs
//
int comp_double(const void *m1, const void *m2)
{
	double x1 = *((double*)m1);
	double x2 = *((double*)m2);

	if (x1 > x2) {
		return -1;
	} else if (x1 == x2) {
		return 0;
	}

	return 1;
};

// main fitting function
//
// sort the datapoints, then iteratively censor them (see above for censoring)
// and re-fit the parameters, until convergence
//
// the parameters are:
//
//      n  ...  number of targets
//     xs  ...  query-target scores
//      p  ...  {lambda, kappa} = paramaters of the distribution
//
// returns the number of unrelated matches, n_unrel
//
long int iterate_censor_fit(long int n, double *xs, double *p)
{
	double old_lambda;
	long   int n_top, n_fit, n_bot, i;

	qsort((void*)xs, n, sizeof(double), comp_double);

	p[0] = old_lambda = 1.0;
	p[1] = 0.0;

	n_top = 0;
	n_fit = (long int)(FR_FIT * ((double)n));
	n_bot = n - n_fit;


	for (i = 0; i < 20; i++) {
		censor(n, xs, &n_top, &n_fit, &n_bot, p);

		conj_grad(n_top, n_fit, n_bot, xs, p);

		if (fabs(old_lambda - p[0]) < 1.0e-10) {
			break;
		}

		old_lambda = p[0];
	};

	return 1 + n_fit + n_bot;
}

// run iterate_censor_fit on a list of matches
//
// this is the main wrapper routine for calls from other parts of PRC
//
void fit(PARAMS *p, MATCH *list)
{
	double *xs;
	MATCH  *match;
	long   int i;

	malloc_1D_array(xs, double, p->n_models);

	match = list;
	for (i = 0; i < p->n_models; i++) {
		xs[i] = (double) match->reverse;
		match = match->next;
	};

	p->n_unrel = iterate_censor_fit(p->n_models, xs, p->p);

	free_unless_null(xs);
}

#ifdef STAND_ALONE
int main(int argc, char **argv)
{
	double *xs, p[2];
	long   int n;

	read_file(argv[1], &n, &xs);
	iterate_censor_fit(n, xs, p);

	printf("lambda=%f\tkappa=%f\n", p[0], p[1]);
}
#endif

