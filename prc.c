/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         PRCx, the profile comparer eXtended, version 1.0.0

		  prc.c: "wrapper" file for the actual binary

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
#include <string.h>
#include <math.h>
#include <time.h>

#include "upgma.h" // Needed for the phylogenitic tree building
#include "ssw.h"   // Needed for fast prefiltering
#include "prc.h"
#include "cs219.h" // The discrete alphabet scoring matrix

// It here is no match or it is negative use this distance for the distance matrix
const double HUGE_DIST = 0.001;

char usage[] =
    "Usage: prc [-options] <hmm> <library> <output>\n"
    "   or  prc [-options] <hmm1> <hmm2>\n"
    "   or  prc [-options] -tree <library> <output>\n"
    "\n"
    "Available options:\n"
    "\n"
    "  -algo  <> : forward (default), viterbi\n"
    "  -MMfn  <> : match-match scoring function; dot1, dot2 (default)\n"
    "  -mode  <> : local-local (default), local-global, global-local,\n"
    "              global-global\n"
    "  -align <> : alignment style; none (default), prc1, prc2, sam1, sam2\n"
    "  -stop  <> : stop looking for more hits when simple<stop (default varies)\n"
    "  -hits  <> : stop looking for more hits when hit_no>hits (default varies)\n"
    "  -Emax  <> : only report hits with E-value<=Emax (default: 10)\n"
    "  -tree     : create a distance matrix and a newick tree of the library\n"
    "              using UPGMA. For a -tree run, hits = 1 and align = none.\n"
    "  -filt     : only check the top 1000 matches as determined using a fast\n"
    "              alignment using Smith-Waterman and a discrete alphabet\n"
    "              only useful if using .prc files (default: off)\n"
    "\n"
    "(See README.txt for more details.)\n"
    "\n"
    "Recognized file extensions:\n"
    "\n"
    KNOWN_EXTENSIONS
    "&\n"
    "  .lib     : library of models (only for second argument!)\n"
    "\n"
    "The library files should be simple text files listing one filename per "
    "line.\n"
    "\n";

// need to have usage[] for common.c
#include "common.c"

char   command_line[1000];
time_t start_time;


// print
//
//   # Command     : ...
//   #    :
//   # Started     : ...
//   #
//
void print_settings(FILE *stream, PARAMS *p, char *init)
{
	struct tm *tm   = localtime(&start_time);
	char *nice_Emax;
	char *aligns[]  = { "none", "prc1", "prc2", "sam1", "sam2" };
	char *MMfn[]    = { "dot1", "dot2" };


	fprintf(stream,

	        "%sCommand     : %s\n"
	        "%sAlgorithm   : %s\n"
	        "%sMatch-match : %s\n"
	        "%sAlign mode  : %s-%s\n"
	        "%sAlignments  : %s\n"
	        "%sSimple stop : %.1f\n"
	        "%sMax hits    : %ld\n",

	        init, command_line,
	        init,
	        (p->algorithm == PARAM_FORW_BACK) ? "forward/backward" : "Viterbi",
	        init, MMfn[p->MM_function],
	        init,
	        (p->align_mode1 == LOCAL) ? "local" : "global",
	        (p->align_mode2 == LOCAL) ? "local" : "global",
	        init, aligns[p->align_style],
	        init, p->stop,
	        init, p->max_hits);

	if (p->E_values) {
		nice_Emax  = nice_string(p->Emax);
		fprintf(stream, "%sMax E-value : %s\n", init, nice_Emax);
		free_unless_null(nice_Emax);
	};

	fprintf(stream,

	        "%sStarted     : %s"
	        "%s\n",

	        init, asctime(tm),
	        init);
}

void print_distrib(FILE *stream, PARAMS *p, char *init)
{
	fprintf(stream,

	        "%sE-value fn  : n_unrel / (1.0 + exp(lambda*reverse + kappa))\n"
	        "%s .. n_unrel : %ld\n"
	        "%s .. lambda  : %f\n"
	        "%s .. kappa   : %f\n",

	        init,
	        init, p->n_unrel,
	        init, p->p[0],
	        init, p->p[1]);

	if (p->n_unrel < 1000)
		fprintf(stream,
		        "%sWarning     : n_unrel < 1000 means poor statistics!\n",
		        init);

	fprintf(stream, "%s\n", init);
}

// print
//
//   # hmm1 start1 end1 ...
//
void print_column_names(FILE *stream, PARAMS *p)
{
	fprintf(stream,
	        "# hmm1\tstart1\tend1\tlength1\thit_no"
	        "\thmm2\tstart2\tend2\tlength2");

	if (p->algorithm == PARAM_VITERBI) {
		fprintf(stream, "\tsimple\treverse");
	} else if (p->algorithm == PARAM_FORW_BACK) {
		fprintf(stream, "\tco-emis\tsimple\treverse");
	};

	if (p->E_values) {
		fprintf(stream, "\t E-value\n");
	} else {
		fprintf(stream, "\n");
	}
}

// left-pad a string to a new length
//
void left_pad(char *str, int new_len)
{
	int i, len, diff;

	len  = strlen(str);
	diff = new_len - len;

	if (diff > 0) {
		for (i = len; i >= 0; i--) {
			str[diff + i] = str[i];
		}

		for (i = diff - 1; i >= 0; i--) {
			str[i] = ' ';
		}
	};
}

// print
//
//   0016282	72	90	124     ...
//
void print_match(FILE *stream, MATCH *m, PARAMS *p)
{
	char *nice_Evalue = nice_string(m->E_value);

	left_pad(nice_Evalue, 8);

	fprintf(stream,
	        "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d",
	        m->hmm1i->name,
	        m->proper.start1 - 1,
	        m->proper.end1 - 1,
	        m->hmm1i->M,
	        m->match_n,
	        m->hmm2i->name,
	        m->proper.start2 - 1,
	        m->proper.end2 - 1,
	        m->hmm2i->M);

	if (p->algorithm == PARAM_FORW_BACK) {
		fprintf(stream, "\t%6.1f", m->sum);
	}

	fprintf(stream,
	        "\t%6.1f\t%6.1f",
	        m->simple,
	        m->reverse);

	if (p->E_values) {
		fprintf(stream, "\t%s\n", nice_Evalue);
	} else {
		fprintf(stream, "\n");
	}

	free_unless_null(nice_Evalue);
}

// print
//
//   >HMM2~hit_no region1=5-100 region2=10-60 reverse=5.0 E-value=0.012
//   alignment_string
//
void print_align(FILE *stream, PARAMS *p, MATCH *m)
{
	char *nice_Evalue;

	if (p->library) {
		// we're saving this to a file -> print the header

		fprintf(stream,
		        ">%s~%s~%d reg1=%d-%d/%d reg2=%d-%d/%d rev=%.1f",
		        m->hmm1i->name,
		        m->hmm2i->name,
		        m->match_n,
		        m->proper.start1 - 1,
		        m->proper.end1 - 1,
		        m->hmm1i->M,
		        m->proper.start2 - 1,
		        m->proper.end2 - 1,
		        m->hmm2i->M,
		        m->reverse);

		if (p->E_values) {
			nice_Evalue = nice_string(m->E_value);
			fprintf(stream, " E=%s", nice_Evalue);
			free_unless_null(nice_Evalue);
		};

		fprintf(stream, "\n");
	} else {
		// this is going to stdout, so a blank to separate from scores
		fprintf(stream, "\n");
	};

	fprintf(stream, "%s", m->alignment);

	if (!p->library) {
		// this is going to stdout, so a blank to separate from scores
		fprintf(stream, "\n");
	};
}

// print out all the matches in p->matches
//
void print_matches(PARAMS *p)
{
	FILE  *scores, *aligns;
	MATCH *match;
	long  int i;

	// set up scores & aligns
	if (p->library) {
		open_file_or_die(scores,  p->output_scores, "w");

		if (p->align_style != ALIGN_NONE) {
			open_file_or_die(aligns,  p->output_aligns, "w");
		} else {
			aligns = stdout;
		};

		print_header(scores, "# ");
		print_settings(scores, p, "# ");

		if (p->E_values) {
			print_distrib(scores, p, "# ");
		}
	} else {
		scores = stdout;
		aligns = stdout;
	};

	print_column_names(scores, p);

	for (i = 0; i < p->n_matches; i++) {
		match = p->matches[i];

		if (SIG_MATCH(p, match)) {
			print_match(scores, match, p);

			if (p->align_style != ALIGN_NONE) {
				print_align(aligns, p, match);
			}
		};
	};

	// close scores & aligns
	fprintf(scores, "# END\n");

	if (p->library) {
		fclose(scores);

		if (p->align_style != ALIGN_NONE) {
			fclose(aligns);
		}

		printf("All done!\n\n");
	};
}

// run two profile files against each other
//
void run_two_files(PARAMS *p)
{
	HMM   *hmm1 = read_HMM(p->model_file1), *hmm2 = read_HMM(p->model_file2);
	HMM   *rev1 = reverse_HMM(hmm1);
	MATCH *matches;

	print_settings(stdout, p, "");

	// 0,1 = "multiple hits, do alignments"
	matches = run_two_HMMs(p, hmm1, rev1, hmm2, 0, 1);

	copy_sort_matches(p, matches);
	print_matches(p);

	free_HMMinfo(hmm1->i);
	free_HMM(hmm1);
	free_HMM(rev1);

	free_HMMinfo(hmm2->i);
	free_HMM(hmm2);
}

// run a profile against a library of profiles
//
// (the library is assumed to be a list of files, one per line)
//
void run_file_against_lib(PARAMS *p)
{
	char  **libfile;
	HMM   *hmm1 = read_HMM(p->model_file1);
	HMM   *rev1 = reverse_HMM(hmm1);
	HMM   *hmm[500];
	long  int i, n, j, progr_ind = 0;
	MATCH *match, *sig_match, *pass1_list = NULL, *pass2_list = NULL;
	size_t n_models = 0;
	int temp;
	int *scores = NULL;
	int *index = NULL;

	/* Read the library file */
	libfile = prc_stringfile(p->model_file2, &n_models);
	if (libfile == NULL) {
		return;
	}

	print_settings(stdout, p, "");

	// prefilter using a fast Smith-Waterman and a discrete alphabet (cs219.lib)
	
	/* If we are prefiltering we need the query to have a discrete alphabet, if 
	 * it does not have one, we need to make it. we then use this to make the query
	 */
	if (p->filt) {
		printf("Filtering   : ");
		/* Needed for prefiltering */
		CONTEXT_LIB *cs219 = NULL;
		s_profile *query = NULL;
		s_align *alignment = NULL;

		/* If we don't have a discrete alphabet yet, make one for the query */
		if (hmm1->discrete_alphabet == NULL) {
			// Load the context library 
			CONTEXT_LIB *cs219 = NULL;
			if ((cs219 = read_context_library("cs219.lib")) == NULL) {
				fprintf(stderr, "Error reading context library from: %s\n", "cs219.lib");
				exit(1);
			}

			// Translate the probabilities to a discrete alphabet
			hmm1->discrete_alphabet = profile_to_discrete_alphabet((double **)&hmm1->Pmat[2], hmm1->i->M, cs219);
			if (hmm1->discrete_alphabet == NULL) {
				fprintf(stderr, "Error converting HMM emissions to discrete alphabet\n");
				exit(1);			
			}
		}
		// Build fast the query profile
		query = ssw_init(hmm1->discrete_alphabet, hmm1->i->M, mat, 219, 2);

		// Allocate space for the scores
		calloc_1D_array(scores, int, n_models);
		calloc_1D_array(index, int, n_models);

		/* Go through all the sequences and calculate their score */
		for (i = 0; i < n_models; i++) {
			if (libfile[i] == NULL) {
				fprintf(stderr, "Please remove the empty line from %s at line %li!\n", p->model_file2, i);
				return;
			}
			// Go through the list 500 at a time
			for (n = 0; (n < 500) && (i < n_models); n++, i++) {
				// Read the HMM
				hmm[n] = read_HMM(libfile[i]);
				if (hmm[n]->discrete_alphabet == NULL) {
					
					if (cs219 == NULL) {
						if ((cs219 = read_context_library("cs219.lib")) == NULL) {
							fprintf(stderr, "Error reading context library from: %s\n", "cs219.lib");
							exit(1);
						}
					}
					hmm[n]->discrete_alphabet = profile_to_discrete_alphabet((double **)&hmm[n]->Pmat[2], hmm[n]->i->M, cs219);
					if (hmm1->discrete_alphabet == NULL) {
						fprintf(stderr, "Error converting HMM emissions to discrete alphabet\n");
						exit(1);			
					}
				}
				// Perform the actual alignment and scoring
				//alignment = ssw_align(query, hmm[n]->discrete_alphabet, hmm[n]->i->M, 90, 30, 0xff, 0xff, 0, 15);
				alignment = ssw_align(query, hmm[n]->discrete_alphabet, hmm[n]->i->M, 90, 30, 0x00, 0x00, 0, 0);

				scores[i] = alignment->score1;
				index[i] = i; //need it later to load the HMM files in the order best->worst
				
				//fprintf(stderr, "%li: %i %s\n", i, alignment->score1, libfile[i]);

				// no longer needed
				align_destroy(alignment);
				free_HMMinfo(hmm[n]->i);
				free_HMM(hmm[n]);

				if (++progr_ind % 100 == 0) {
					progr_ind = 0;
					printf(".");
					fflush(stdout);
				};
			}
		}

		if (cs219) {
			free_context_library(cs219);
		}
		init_destroy(query);
		
		// Sort the scores from hight to low
		for(i = 0; i < n_models; i++) {
			for(j = i + 1; j < n_models; j++) {
				if(scores[i] < scores[j]) {
					// sort the scores
					temp = scores[i];
					scores[i] = scores[j];
					scores[j] = temp;
					// and the index of the scores
					temp = index[i];
					index[i] = index[j];
					index[j] = temp;
				}
			}
		}

		//int pass = 1010; // Only allow enough sequences to pass to calculate a accurate E-value
		int pass = 1010; // Only allow enough sequences to pass to calculate a accurate E-value
		if (n_models > pass) {
			n_models = pass;
		}

		// printf("Scores:\n");
		// for(i = 0; i < n_models; i++) {
		// 	printf("%i: %i %s\n", index[i], scores[i], libfile[index[i]]);
		// }
		// printf("\n");

		printf(" Done!\n\n");
	}

	// pass 1: - if E-values, get the top match for each model in the library
	//         - if no E-values, just get all the matches straight away
	//
	printf("Scoring     : ");
	fflush(stdout);

	for (i = 0; i < n_models; i++) {
		// burst read of HMM files to maximize caching
		for (n = 0; (n < 500) && (i < n_models); n++, i++) {
			if (p->filt) {
				hmm[n] = read_HMM(libfile[index[i]]);	
				p->n_models++; // needed here for statistics do NOT move!
			} else {
				hmm[n] = read_HMM(libfile[i]);
				p->n_models++; // needed here for statistics do NOT move!
			}
		};

		// and now take your time running them
		for (j = 0; j < n; j++) {
			if (p->E_values) {
				if (p->max_hits > 1) {
					// 1,0 = "single hit, no alignments"
					match = run_two_HMMs(p, hmm1, rev1, hmm[j], 1, 0);
				} else {
					// 1,1 = "single hit, alignments"
					match = run_two_HMMs(p, hmm1, rev1, hmm[j], 1, 1);
				};
			} else {
				// 0,1 = "multiple hits, alignments"
				match = run_two_HMMs(p, hmm1, rev1, hmm[j], 0, 1);
			};

			if (match == NULL) {
				// if no match, get rid of the info now
				free_HMMinfo(hmm[j]->i);
			} else {
				// otherwise, tie the info to the first match
				match->free_hmm2i = 1;
			};

			// NB this doesn't affect the info
			free_HMM(hmm[j]);

			// keep all the matches in a list starting at pass1_list
			pass1_list = join_lists(match, pass1_list);

			if (++progr_ind % 100 == 0) {
				progr_ind = 0;
				printf(".");
				fflush(stdout);
			};
		};
	};

	printf(" Done!\n\n");
	fflush(stdout);

	// sanity check
	if (p->E_values && (p->n_models < 10)) {
		printf("Error: Ahem, you have %ld model(s) in your library. "
		       "Can only do statistics on 10 models or more!\n\n",
		       p->n_models);

		p->E_values = 0;
	};

	if (p->E_values) {
		// fit the distribution to the matches in pass1_list
		fit(p, pass1_list);
		print_distrib(stdout, p, "");

		for (match = pass1_list; match; match = match->next) {
			match->E_value = E_value(p->n_unrel, match->reverse, p->p);
		}
	};

	// pass 2: get all hits for models that make good matches
	if ((p->E_values) && (p->max_hits > 1)) {
		printf("Re-scoring  : ");
		fflush(stdout);
		progr_ind = 0;

		for (match = pass1_list; match;) {
			// burst read of 500 models to maximize caching
			for (i = 0; (i < 500) && match; match = match->next)
				if (WEAK_MATCH(p, match)) {
					hmm[i++] = read_HMM(match->hmm2i->filename);
				}

			// now run the 500 models
			for (j = 0; j < i; j++) {
				// 0,1 = "multiple hits, alignments"
				sig_match = run_two_HMMs(p, hmm1, rev1, hmm[j], 0, 1);

				if (sig_match) {
					// tie hmm2i to the first match
					sig_match->free_hmm2i = 1;
				} else {
					// get rid of hmm2i now, as it's no longer needed
					free_HMMinfo(hmm[j]->i);
				};

				free_HMM(hmm[j]);
				pass2_list = join_lists(sig_match, pass2_list);

				if (++progr_ind % 100 == 0) {
					progr_ind = 0;
					printf(".");
					fflush(stdout);
				};
			};
		};
		printf(" Done!\n\n");
		fflush(stdout);
			
		copy_sort_matches(p, pass2_list);
		print_matches(p);
		free_matches(pass2_list);
	} else {
		// no 2nd pass, print the results from the first pass
		copy_sort_matches(p, pass1_list);
		print_matches(p);
	};

	free(libfile);
	free_matches(pass1_list);
	free_HMMinfo(hmm1->i);
	free_HMM(hmm1);
	free_HMM(rev1);
}

// run a library of profiles against itself and construct a distance
// matrix and a tree using UPGMA
//
// (the library is assumed to be a list of files, one per line)
//
void run_lib(PARAMS *p)
{
	char **libfile = NULL;
	int progr_ind = 0;
	HMM   *hmm1 = NULL;
	//HMM   *rev1 = NULL;
	HMM   *hmm2 = NULL;
	double distance = 0.0;
	long int i, n;
	size_t n_models = 0;
	MATCH *match;
	matrix_t *distmat;
	tree_t tree;
	char **labels;

	/* Read the library file */
	libfile = prc_stringfile(p->model_file1, &n_models);
	if (libfile == NULL) {
		return;
	}

	if (n_models < 3) {
		die1("It is impossible to make a tree with less than 3 sequences!\n");
	}

	/* Allocate space for the names of the models */
	labels = (char **) calloc(n_models, sizeof(char**));

	/* Allocate a distance matrix */
	prc_newmat(&distmat, n_models, n_models);

	print_settings(stdout, p, "");

	printf("Scoring     : ");
	fflush(stdout);

	/* Create the distance matrix */
	for (i = 0; i < n_models; i++) {
		if (libfile[i] == NULL) {
			fprintf(stderr, "Please remove the empty line from %s at line %li!\n", p->model_file1, i);
			return;
		} 
		hmm1 = read_HMM(libfile[i]);
		// Actually we do not use the reverse. Maybe we should remove it? Or use it somehow?
		// For now we removed it. Small tests indicate that without the additional
		// reverse computation, a comparison of 248 pHHMS takes 23 min versus 36 min.
		// This difference is obviously less for smaller comparisons, but more for larger ones!
		//rev1 = reverse_HMM(hmm1);
		distmat->data[i][i] = 0.0;
		/* Copy the name of the HMM */
		labels[i] = strdup(hmm1->i->name);

		for (n = i + 1; n < n_models; n++) {
			if (libfile[n] == NULL) {
				fprintf(stderr, "Please remove the empty line from %s at line %li!\n", p->model_file1, n);
				return;
			} 
			hmm2 = read_HMM(libfile[n]);
			// We are not doing anything useful with the reverse scores
			//match = run_two_HMMs(p, hmm1, rev1, hmm2, 1, 0);
			match = simple_run_two_HMMs(p, hmm1, hmm2, 1, 0);

			if (match) {
				/* Save the match in the distmat */
				if (match->simple > 0.0) {
					distance = 1.0 / match->simple;
					distmat->data[i][n] = distmat->data[n][i] = distance;
				} else {
					distance = 1.0 / HUGE_DIST;
					distmat->data[i][n] = distmat->data[n][i] = distance;
				}
				/* We no longer need the match */
				free_matches(match);
			} else {
				/* if there is no match we make it so there is a "huge" distance
				 * distance = 1.0 / 0.001; // HUGE_DIST defined at the top of this file
				 */
				distance = 1.0 / HUGE_DIST;
				distmat->data[i][n] = distmat->data[n][i] = distance;
			}

			free_HMMinfo(hmm2->i);
			// NB this doesn't affect the info
			free_HMM(hmm2);

			// Show that we are making progress
			if (++progr_ind % 100 == 0) {
				progr_ind = 0;
				printf(".");
				fflush(stdout);
			}
		}
		free_HMMinfo(hmm1->i);
		free_HMM(hmm1);
		//free_HMM(rev1);
	}
	
	printf(" Done!\n\n");
	fflush(stdout);

	free(libfile);

	
	printf("Saving Distance matrix...");
	fflush(stdout);
	FILE *outdist;
	open_file_or_die(outdist, p->output_distmat, "w");
	save_distmat(outdist, distmat, labels);
	fclose(outdist);

	printf(" Done!\n\n");

	/* Create the UPGMA tree */
	printf("Creating Tree...");
	fflush(stdout);	
	
	/* Perform the UPGMA clustering using average linkage */
	prc_upgma(&tree, distmat, LINKAGE_AVG, labels);

	printf(" Done!\n\n");

	printf("Saving Tree...");
	fflush(stdout);
	
	/* Save the tree */
	FILE *outtree;
	open_file_or_die(outtree, p->output_tree, "w");
	save_tree(outtree, &tree);
	fclose(outtree);
	
	printf(" Done!\n\n");
	fflush(stdout);
	
	free_tree(&tree);
	prc_freemat(distmat);

	/* Free the labels (HMM names) */
	for (i = 0; i < n_models; i++) {
		free(labels[i]);
	}
	free(labels);
}

// does the string say "local" or "global", and similar questions
//
#define TYPE_MODE   0
#define TYPE_ALGO   1
#define TYPE_MMFN   2
#define TYPE_ALIGN  3
int parse_param(int type, char *str)
{
	char *error_str[4] = { "alignment mode",
	                       "algorithm",
	                       "scoring function",
	                       "alignment style"
	                     };

	int  number[4]     = { 2, 2, 2, 5 };

	char *param[4][5]  = { { "local", "global", "", "", "" },
	                       { "forward", "viterbi", "", "", "" },
	                       { "dot1", "dot2", "", "", "" },
	                       { "none", "prc1", "prc2", "sam1", "sam2" }
	};

	int  length[4][5]  = { { 5, 6, 0, 0, 0 },
	                       { 7, 7, 0, 0, 0 },
	                       { 4, 4, 0, 0, 0 },
	                       { 4, 4, 4, 4, 4 }
	};

	int  ret_val[4][5] = { { LOCAL, GLOBAL, -1, -1, -1 },
		{ PARAM_FORW_BACK, PARAM_VITERBI, -1, -1, -1 },
		{ DOT1, DOT2, -1, -1, -1 },
		{ ALIGN_NONE, ALIGN_PRC1, ALIGN_PRC2, ALIGN_SAM1, ALIGN_SAM2 }
	};

	int  i;


	if ((type < 0) || (type > 4)) {
		die2("Uknown parameter type #%d", type);
	}

	if (!str) {
		arg_error("Error: Uknown %s ''!", error_str[type]);
	}

	for (i = 0; i < number[type]; i++) {
		if (strncmp(str, param[type][i], length[type][i] + 1) == 0) {
			return ret_val[type][i];
		}
	};

	arg_error("Error: Uknown %s '%s'!", error_str[type], str);

	// to keep GCC ehappy
	return 0;
}

// nowhere else to put it, so stuck these two routines here
//
PARAMS* alloc_params(void)
{
	PARAMS *p;

	calloc_1D_array(p, PARAMS, 1);

	// -1 = not set yet
	// (for stop this is handled via stop_flag in main())
	p->algorithm   = -1;
	p->MM_function = -1;
	p->align_mode1 = -1;
	p->align_mode2 = -1;
	p->align_style = -1;
	p->E_values    = -1;
	p->Emax        = -1.0;
	p->max_hits    = -1;
	p->tree        =  0; // build a distance matrix and tree? 0 = no / 1 = yes
	p->filt        =  0; // use fast smith-waterman prefilter? 0 = no / 1 = yes

	return p;
}

// fully free a PARAMS structure
//
// also frees all the matches in p->matches
//
void free_params(PARAMS *p)
{
	free_unless_null(p->matches);
	free_unless_null(p->output_scores);
	free_unless_null(p->output_aligns);
	free_unless_null(p->model_file1);
	free_unless_null(p->model_file2);
	free_unless_null(p);
}

// the main routine: parse the command-line parameters
//
int main(int argc, char **argv)
{
	PARAMS *p;
	int    argi, len, stop_flag = -1;
	FILE   *check_file;


	// record the start time
	start_time = time(NULL);

	// sort out the version string
	sprintf(version, "PRCx " VERSION " (%s, %s, %s), compiled on %s",
#if   PROF_HMM_TRANS==PLAN9
	        "PLAN9"
#elif PROF_HMM_TRANS==PLAN7
	        "PLAN7"
#endif
	        ,
#if   PAIR_HMM_STATES==SPACE9
	        "SPACE9"
#elif PAIR_HMM_STATES==SPACE8
	        "SPACE8"
#elif PAIR_HMM_STATES==SPACE6
	        "SPACE6"
#elif PAIR_HMM_STATES==SPACE5
	        "SPACE5"
#endif
	        ,
#if   PAIR_HMM_TRANS==ALL_TRANS
	        "ALL_TRANS"
#elif PAIR_HMM_TRANS==VIA_MM
	        "VIA_MM"
#endif
	        , __DATE__);

	// print out the header
	print_header(stdout, "");

	// keep the command line
	command_line[0] = 0;
	for (argi = 0; argi < argc; argi++) {
		strcat(command_line, argv[argi]);
		strcat(command_line, " ");
	};

	// sort out the options
	p = alloc_params();

	for (argi = 1; argi < argc; argi++) {
		// onto the filenames ...
		if (argv[argi][0] != '-') {
			break;
		}

		if (strncmp(argv[argi], "-algo", 6) == 0) {
			if (p->algorithm != -1) {
				arg_error("Error: Alignment algorithms set more than once!");
			}

			p->algorithm = parse_param(TYPE_ALGO, argv[++argi]);
		} else if (strncmp(argv[argi], "-MMfn", 6) == 0) {
			if (p->MM_function != -1) {
				arg_error("Error: Scoring function set more than once!");
			}

			p->MM_function = parse_param(TYPE_MMFN, argv[++argi]);
		} else if (strncmp(argv[argi], "-mode", 6) == 0) {
			if ((p->align_mode1 != -1) || (p->align_mode2 != -1)) {
				arg_error("Error: Alignment mode set more than once!");
			}

			p->align_mode1 = parse_param(TYPE_MODE, strtok(argv[++argi], "-"));
			p->align_mode2 = parse_param(TYPE_MODE, strtok(NULL,         "-"));
		} else if (strncmp(argv[argi], "-align", 8) == 0) {
			if (p->align_style != -1) {
				arg_error("Error: Alignment style set more than once!");
			}

			p->align_style = parse_param(TYPE_ALIGN, argv[++argi]);
		} else if (strncmp(argv[argi], "-Emax", 6) == 0) {
			if (p->Emax != -1.0) {
				arg_error("Error: Emax threshold set more than once!");
			}

			if (sscanf(argv[++argi], "%lf", &p->Emax) != 1) {
				arg_error("Error parsing the Emax cutoff '%s'!", argv[argi]);
			}

			if (p->Emax < 0.0)
				arg_error("Error: Emax cutoff '%s' seems to be negative!",
				          argv[argi]);
		} else if (strncmp(argv[argi], "-stop", 6) == 0) {
			if (stop_flag != -1) {
				arg_error("Error: Stop threshold set more than once!");
			}

			if (sscanf(argv[++argi], "%lf", &p->stop) != 1) {
				arg_error("Error parsing the stop threshold '%s'!", argv[argi]);
			}

			stop_flag = 1;
		} else if (strncmp(argv[argi], "-hits", 5) == 0) {
			if (p->max_hits != -1) {
				arg_error("Error: Maximum number of hits set more than once!");
			}

			if (sscanf(argv[++argi], "%ld", &p->max_hits) != 1) {
				arg_error("Error parsing max number of hits '%s'!", argv[argi]);
			}

			if (p->max_hits < 1) {
				arg_error("Error: The maximum number of hits must be > 0 !");
			}
		} else if ((strncmp(argv[argi], "-h", 3) == 0) ||
		           (strncmp(argv[argi], "--help", 7) == 0)) {
			print_help_die();
		} else if (strncmp(argv[argi], "-tree", 5) == 0) {
			if (p->tree != 0) {
				arg_error("Error: Tree option set more than once!");
			}

			p->tree = 1;
		} else if (strncmp(argv[argi], "-filt", 5) == 0) {
			if (p->filt != 0) {
				arg_error("Error: Filtering option set more than once!");
			}

			p->filt = 1;
		} else {
			arg_error("Error: Unknown option '%s'!", argv[argi]);
		};
	};

	// set the defaults for anything that wasn't set
	if (p->algorithm   == -1) {
		p->algorithm   = PARAM_FORW_BACK;
	}
	if (p->MM_function == -1) {
		p->MM_function = DOT2;
	}
	if (p->align_mode1 == -1) {
		p->align_mode1 = LOCAL;
	}
	if (p->align_mode2 == -1) {
		p->align_mode2 = LOCAL;
	}
	if (p->align_style == -1) {
		p->align_style = ALIGN_NONE;
	}
	if (p->Emax        == -1.0) {
		p->Emax        = 10.0;
	}

	// set p->stop & p->max_hits
	if ((p->align_mode1 == LOCAL) && (p->align_mode2 == LOCAL)) {
		if (stop_flag == -1) {
			p->stop = 1.5;
		}

		// can't set p->max_hits because don't know yet whether this run is a
		// pairwise or library one
	} else if (((p->align_mode1 == LOCAL) && (p->align_mode2 == GLOBAL)) ||
	           ((p->align_mode1 == GLOBAL) && (p->align_mode2 == LOCAL))) {
		if (p->algorithm == PARAM_FORW_BACK) {
			if (stop_flag == -1) {
				p->stop = 1.5;
			}
		} else if (p->algorithm == PARAM_VITERBI) {
			if (stop_flag == -1) {
				p->stop = 1.5;
			}
		} else {
			die2("Unknown algorithm #%d!", p->algorithm);
		};

		if (p->max_hits == -1) {
			p->max_hits = 10;
		}
	} else if ((p->align_mode1 == GLOBAL) && (p->align_mode2 == GLOBAL)) {
		if (p->algorithm == PARAM_FORW_BACK) {
			if (stop_flag == -1) {
				p->stop = -2.0;
			}
		} else if (p->algorithm == PARAM_VITERBI) {
			if (stop_flag == -1) {
				p->stop = -2.0;
			}
		} else {
			die2("Unknown algorithm #%d!", p->algorithm);
		};
	} else {
		die3("Uknown alignment mode %d-%d!", p->align_mode1, p->align_mode2);
	};


	// and now the filenames
	if (argc - argi == 0) {
		print_help_die();
	} else if (((argc - argi == 1) && !p->tree ) || (argc - argi > 3)) {
		arg_error("Error: Incorrect number of arguments!");
	}

	// Tree building is the only mode that only needs one file, the library and pairwise need two
	if (p->tree) {
		len = strlen(argv[argi]);
		calloc_1D_array(p->model_file1, char, len + 1);
		strncpy(p->model_file1, argv[argi], len + 1);

		if (strncmp(&p->model_file1[len - 4], ".lib", 4) != 0) {
			arg_error("Error: Library (.lib) is needed for tree building not '%s'!", p->model_file1);
		} else {
			p->library = 1;
		}

		/* We only use the best hit for the distance matrix */
		if (p->max_hits == -1) {
			p->max_hits = 1;
		} else if (p->max_hits > 1) {
			arg_error("Error: For tree building we only use 1 hit, max hits of '%i' ignored !", p->max_hits);
			p->max_hits = 1;
		}

		/* We are not going to calculate an E-value */
		/* We only use the best hit for the distance matrix */
		if (p->E_values == -1) {
			p->E_values = 0;
		} else if (p->E_values > -1) {
			arg_error("Error: For tree building we only use the simple score, E-value argument ignored !");
			p->E_values = 0;
		}

		/* We need atleast the library file and the output */
		if (argc - argi != 2) {
			arg_error("Error: Incorrect number of arguments!");
		}

		// output
		len = strlen(argv[++argi]);
		calloc_1D_array(p->output_distmat, char, len + 9);
		calloc_1D_array(p->output_tree, char, len + 8);
		sprintf(p->output_distmat, "%s.distmat", argv[argi]);
		sprintf(p->output_tree, "%s.newick", argv[argi]);

		// make sure that the files don't already exist
		check_file = fopen(p->output_distmat, "r");
		if (check_file != NULL) {
			fclose(check_file);
			arg_error("Error: The file '%s' already exists. "
			          "Please remove it first!", p->output_distmat);
		}

		check_file = fopen(p->output_tree, "r");
		if (check_file != NULL) {
			fclose(check_file);
			arg_error("Error: The file '%s' already exists. "
			          "Please remove it first!", p->output_tree);
		}
		run_lib(p);
	} else {
		// first file
		len = strlen(argv[argi]);
		calloc_1D_array(p->model_file1, char, len + 1);
		strncpy(p->model_file1, argv[argi], len + 1);
		if (!known_extension(p->model_file1)) {
			arg_error("Error: Unknown extension for model file '%s'!", p->model_file1);
		}

		// second file or library
		len = strlen(argv[++argi]);
		calloc_1D_array(p->model_file2, char, len + 1);
		strncpy(p->model_file2, argv[argi], len + 1);

		if (strncmp(&p->model_file2[len - 4], ".lib", 4) == 0) {
			// the second argument is a library of HMMs
			p->library = 1;
			if ((p->align_mode1 == LOCAL) && (p->align_mode2 == LOCAL)) {
				p->E_values = 1;
			} else {
				p->E_values = 0;
			}

			// p->max_hits only unset here for local-local runs
			if (p->max_hits == -1) {
				if (p->E_values) {
					p->max_hits = 100000;
				} else {
					p->max_hits = 100;
				}
			};

			if (argc - argi != 2) {
				arg_error("Error: Incorrect number of arguments!");
			}

			// output
			len = strlen(argv[++argi]);
			calloc_1D_array(p->output_scores, char, len + 8);
			calloc_1D_array(p->output_aligns, char, len + 8);
			sprintf(p->output_scores, "%s.scores", argv[argi]);
			sprintf(p->output_aligns, "%s.aligns", argv[argi]);

			// make sure that the files don't already exist
			check_file = fopen(p->output_scores, "r");
			if (check_file != NULL) {
				fclose(check_file);
				arg_error("Error: The file '%s' already exists. "
				          "Please remove it first!", p->output_scores);
			};

			if (p->align_style != ALIGN_NONE) {
				check_file = fopen(p->output_aligns, "r");
				if (check_file != NULL) {
					fclose(check_file);
					arg_error("Error: The file '%s' already exists. "
					          "Please remove it first!", p->output_aligns);
				};
			};


			run_file_against_lib(p);
		} else if (known_extension(p->model_file2)) {
			// the second argument is another HMM
			p->library  = 0;
			p->E_values = 0;

			// p->max_hits only unset here for local-local runs
			if (p->max_hits == -1) {
				p->max_hits = 100;
			}

			if (argc - argi > 1) {
				arg_error("Error: Incorrect number of arguments!");
			}

			run_two_files(p);
		} else {
			arg_error("Error: Uknown extension for model file '%s'!",
			          p->model_file2);
		};
	}
	free_params(p);

	return 0;
}
