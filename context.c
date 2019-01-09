#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <stdbool.h>

#include "prc.h"

const char *ALPHABET_PRC = "ACDEFGHIKLMNPQRSTVWY";
const char *GAPS = "-.";
const char *ALPHABET_PRC_INC_GAPS = "ACDEFGHIKLMNPQRSTVWY-.X";
const char *ANY = "X";
const char *ANY_OR_GAP = "-.X";

#ifndef MAX
	#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
	#define MIN(x,y)        ((x) < (y) ? (x) : (y))
#endif

#ifdef DBL_DECIMAL_DIG
	#define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else  
	#ifdef DECIMAL_DIG
		#define OP_DBL_Digs (DECIMAL_DIG)
	#else  
		#define OP_DBL_Digs (DBL_DIG + 3)
	#endif
#endif

#define IS_ODD(n) (n & 1)
#define IS_EVEN(n) !IS_ODD(n)

#include <stdbool.h>

/* Extended exponential function */
double eexp(double xin)
{
	if (isnan(xin)) {
		return 0.0;
	} else {
		return exp(xin);
	}
}

/* Extended logarithm function */
double eln(double xin)
{
	if (fabs(xin) < DBL_EPSILON) {
		return NAN;
	} else if (xin > 0) {
		return log(xin);
	} else {
		fprintf(stderr, "Negative input: %lf\n", xin);
		return xin;
	}
}

/* Extended logarithm sum function */
double elnsum(double xin, double yin)
{
	if (isnan(xin) || isnan(yin)) {
		if (isnan(xin)) {
			return yin;
		} else {
			return xin;
		}
	} else if (xin > yin) {
		return xin + eln(1 + exp(yin-xin));
	} else {
		return yin + eln(1 + exp(xin-yin));
	}
}

/* extended logarithm product function */
double elnproduct(double xin, double yin)
{
	if (isnan(xin) || isnan(yin)) {
		return NAN;
	} else {
		return (xin + yin);
	}
}

/* Convert a character belonging to a alphabet to an int */
/* returns -1 if c is not an alphabetic character */
int prc_char2int(char c, const char *alphabet)
{
	char *p = NULL;
	p = strchr(alphabet, toupper(c));

	return p ? (int)(p - alphabet) : -1;
}

/* returns ? if c is not an alphabetic character */
char prc_int2char(int c, const char *const alphabet)
{ 
  if (c < strlen(alphabet)) {
    return alphabet[c]; 
  }

  return '?';
}

inline static double prc_sqr(double xin) {
  return xin*xin;
}

/* Reads an entire file into an array of strings, needs only a single call to free */
char **prc_stringfile2(char *filename, size_t *number_of_lines)
{
    size_t count = 1;
    char **sfile = NULL;
    char *p;

    FILE *f = NULL;
    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Unable to open: %s\n", filename);
        return NULL;
    }

    /* Determine the file size */
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    /* Read the file into a temporary buffer */
    char *buffer = malloc(fsize + 1);
    fread(buffer, fsize, 1, f);
    if (!buffer) {
        fprintf(stderr, "Failed to read %s into memory\n", filename);
        return NULL;
    }
    buffer[fsize] = 0;
    
    /* Close the file */
    fclose(f);

    /* Count the number of new lines */
    p = buffer;
    size_t i = 0;
    while (p[i]) {
        if (p[i] == '\n') {
            if ( p[i+1] == '\r') {
                count++;
                i++;
            } else {
                count++;
            }
        } else if (*p == '\r') {
            count++;
        }
        i++;
    }

    if (number_of_lines) {
        *number_of_lines = count;
    }

    /* Allocate space to keep the entire file */
    sfile = (char **) malloc(sizeof(char *) * (count + 1) + fsize + 1);
    if (!sfile) {
        fprintf(stderr, "Could not copy the data\n");
        return NULL;
    }
    sfile[count] = NULL;
    /* Copy in the original data */
    memcpy(&sfile[count + 1], buffer, fsize + 1);
    
    free(buffer);
    buffer = (char *) &sfile[count + 1];

    /* Go over everything again and set the pointers */
    p = buffer;
    i = 0;
    count = 0;
    sfile[count] = &p[i];
    while (p[i]) {
        if (p[i] == '\n') {
            if ( p[i+1] == '\r') {
                p[i] = '\0';
                p[i+1] = '\0';
                count++;
                i++;
                if (p[i+1]) {
                    sfile[count] = &p[i+1];
                }
            } else {
                p[i] = '\0';
                count++;
                if (p[i+1]) {
                    sfile[count] = &p[i+1];
                }
            }
        } else if (*p == '\r') {
            p[i] = '\0';
            count++;
            if (p[i+1]) {
                sfile[count] = &p[i+1];
            }
        }
        i++;
    }

    return sfile;
}

/* Compute (numerical stable) sum of the data using the Neumaier summation algorithm */
double prc_sum(double *data, int n)
{
	double c = 0.0; /* A running compensation for lost low-order bits. */

	double sum = data[0];
	for (int i = 1; i < n; i++) {
		double t = sum + data[i];
		if (fabs(sum) >= fabs(data[i])) {
			/* If sum is bigger, low-order digits of input[i] are lost. */
			c += (sum - t) + data[i];
		} else {
			/* Else low-order digits of sum are lost */
			c += (data[i] - t) + sum;
		}
		sum = t;
	}

	return sum + c;
}

/* Allocates arbitrary dimensional arrays with one call to calloc.
 * rank = the number of dimensions (e.g. 3 dimensions)
 * shape = an rank sized array indicating the content of each dimension 
 *      (e.g. x =3 y = 9, z = 2 shape[0] = 3 shape[1] = 9 shape[2] = 2)
 * item_size = the size of each element (e.g. sizeof(int))
 */ 
void *prc_allocarray(int rank,int *shape,size_t item_size)
{
  size_t size,i,j,dataSpace,pointerSpace,pointers,nextLevelIncrement;
  char*memory,*pc,*nextpc;

  if (rank < 2) {
    if (rank < 0) {
      fprintf(stderr, "invalid negative rank argument passed to allocarray\n");
      exit(1);
    }
    size = rank < 1 ? 1 : *shape;
    return calloc(size,item_size);
  }

  pointerSpace = 0, dataSpace = 1;
  for (i = 0; i < rank-1; ++i)
    pointerSpace += (dataSpace *= shape[i]);

  pointerSpace *= sizeof(char*);
  dataSpace *= shape[i]*item_size;
  memory = pc = calloc(1,pointerSpace+dataSpace);
  pointers = 1;
  for (i = 0; i < rank-2; ) {
    nextpc = pc + (pointers *= shape[i])*sizeof(char*);
    nextLevelIncrement = shape[++i]*sizeof(char*);
    for (j = 0; j < pointers; ++j)
      *((char**)pc) = nextpc, pc+=sizeof(char*), nextpc += nextLevelIncrement;
  }

  nextpc = pc + (pointers *= shape[i])*sizeof(char*);
  nextLevelIncrement = shape[++i]*item_size;
  for (j = 0; j < pointers; ++j)
    *((char**)pc) = nextpc, pc+=sizeof(char*), nextpc += nextLevelIncrement;
  
  return memory;
}

/* Allocate space for an context library 
 * size = the number of context profiles in the library
 * wlen = the size of the window
 * nalph = the number of characters in the alphabet (= 20 AAs)
 */
CONTEXT_LIB *new_context_library(int size, int wlen, int nalph)
{
	int rank = 3;
	int shape[3];

	CONTEXT_LIB *context_lib = calloc(1, sizeof(CONTEXT_LIB));
	
	/* Default options */
	context_lib->pc_admix         = 0.30; //pca
	context_lib->pc_ali           = 4.0;  //pcb
	context_lib->weight_center    = 1.6;
	context_lib->weight_decay     = 0.85;
	context_lib->weight_as        = 1000.0;
	context_lib->wlen             = wlen;
	context_lib->size             = size;
	context_lib->nalph            = nalph;

	/* We need an odd number to determine the center */
	if (IS_EVEN(context_lib->wlen)) {
		fprintf(stderr, "The window length indicated int the context library should be odd, but is %i\n", context_lib->wlen);
		return NULL;
	}
	/* Calculate the center */
	context_lib->center = (int) ((context_lib->wlen - 1.0) / 2.0);

	context_lib->weights = calloc(wlen, sizeof(double));
	
	context_lib->prior = calloc(size, sizeof(double));

	shape[0] = size;
	shape[1] = wlen;
	shape[2] = nalph;	
	context_lib->probs = prc_allocarray(rank, shape, sizeof(double));

	return context_lib;
}

/* free the conext library struct */
void free_context_library(CONTEXT_LIB *context_lib)
{
	free(context_lib->prior);
	free(context_lib->probs);
	free(context_lib);
}

/* read a context library from a file*/
CONTEXT_LIB *read_context_library(char *filename) 
{
	CONTEXT_LIB *context_lib = NULL;
	char **library = NULL;
	size_t length = 0;
	int i, n;
	int size = 0;
	int wlen = 0;

	if ((library = prc_stringfile2(filename, &length)) == NULL) {
		fprintf(stderr, "Error opening: %s\n", filename);
		return NULL;
	}
	/* Parse the header of the Context Library */
	/* The header should start with this line */
	if (strncmp(library[0], "ContextLibrary", 14)) {
		fprintf(stderr, "%s is not a valid context library\n", filename);
		return NULL;
	}

	/* The next line indicates the number of context profiles in the context library */
	if (strncmp(library[1], "SIZE", 4)) {
		fprintf(stderr, "Unable to parse the SIZE line of the context library\n");
		return NULL;
	} else {
		if (sscanf(&library[1][4],"%i", &size) != 1) {
			fprintf(stderr, "Unable to parse the SIZE of the context library\n");
			return NULL;
		}
	}

	/* The next line indicates the window size of the context profiles in the context library */
	if (strncmp(library[2], "LENG", 4)) {
		fprintf(stderr, "Unable to parse the LENGTH line of the context library\n");
		return NULL;
	} else {
		if (sscanf(&library[2][4],"%i", &wlen) != 1) {
			fprintf(stderr, "Unable to parse the LENGTH of the context library\n");
			return NULL;
		}
	}

	/* Now we have enough information to allocate a new context library structure */
	context_lib = new_context_library(size, wlen, 20);
	context_lib->size = size;
	context_lib->wlen = wlen;
	context_lib->nalph = 20;

	/* The remaining lines contain the contex profiles. Currently, we only parse a 
	 * small subselection of the original hhsuite context profiles.
	 */
	n = 0;
	char *s;
	int score, curr_wlen, curr_aa;
	for (i = 0; i < length; i++) {
		/* Skip empty lines */
		if (library[i] == NULL) {
			fprintf(stderr, "Malformed line in %s at line %i!\n", filename, i);
			continue;
		}

		/* Find the line containing the prior probability of the context profile occuring.
		 * This is stored as a raw probability score. For furthur calculations we convert
		 * this score using the natural logarithm (eln; mainly for numerical stability)
		 */
		if (!strncmp(library[i], "PRIOR", 5)) {
			if (sscanf(&library[i][5],"%lf", &(context_lib->prior[n])) != 1) {
				fprintf(stderr, "Unable to parse the PRIOR at line %i\n", i + 1);
			}
			context_lib->prior[n] = eln(context_lib->prior[n]);
		}

		/* After we got the prior probability, search for the line containing the amino acid
		 * probabilities. These are stored using scaled (by 1000) log2 values as in the original 
		 * hhsuite implementation. So we first convert these values to the raw probabilities and 
		 * then to the natural logarithm (eln). This is again for numerical stability.
		 */
		if (!strncmp(library[i], "PROBS", 5)) {
			/* Move to the next line (the line containing the actual numbers) */
			i++;
			/* Skip the number indicating the position in the window */
			s = strtok(library[i], " \t\n");
			curr_wlen = 0;
			curr_aa = 0;
			/* Every context profile is ended with a // */
			while (library[i][0] != '/') {
				/* get the string containing the first probability */
				s = strtok(NULL, " \t\n");
				/* Convert to an int */
				score = atoi(s);
				/* Convert the scale log2 probability to the natural logarithm */
				context_lib->probs[n][curr_wlen][curr_aa] = eln(pow(2, -(double) score/1000.0));
				/* Move to the next AA probability */
				curr_aa++;
				/* If we read 20 amino acids, we need to move to the next line of the window,
				 * again containing 20 amino acids (so we reset the count)
				 */
				if (curr_aa == 20) {
					curr_aa = 0;
					curr_wlen++;
					/* Move to the next line */
					i++;
					/* Skip the number indicating the position in the window */
					s = strtok(library[i], " \t\n");
				}
			}
			/* We are going to read the next context profile */
			n++;
		}
	}

	/* Done */
	free(library);

	return context_lib;
}

/* Save a context library to a file */
void write_context_library(char *filename, CONTEXT_LIB *context_lib) {
    int i, n = 0, j;
    FILE *fout = NULL;
    fout = fopen(filename, "w");
    if (!fout) {
        fprintf(stderr, "Unable to open: %s\n", filename);
        return;
    }

    // Print header section
    fputs("ContextLibrary\n", fout);
    fprintf(fout, "SIZE\t%i\n", context_lib->size);
    fprintf(fout, "LENG\t%i\n", context_lib->wlen);
    for (i = 0; i < context_lib->size; i++) {
    	fputs("ContextProfile\n", fout);
    	fprintf(fout, "NAME\t%i\n", i + 1);
    	fprintf(fout, "PRIOR\t%.*e\n", OP_DBL_Digs - 1, eexp(context_lib->prior[i]));
    	fputs("PROBS\t", fout);
    	for (j = 0; j < context_lib->nalph; j++) {
    		fprintf(fout, "\t%c", ALPHABET_PRC[j]);
    	}
    	fputs("\n", fout);

    	for (n = 0; n < context_lib->wlen; n++) {
    		fprintf(fout, "%i\t", n + 1);
    		for (j = 0; j < context_lib->nalph; j++) {
    			fprintf(fout, "\t%i", (int) -(1000.0 * log2(eexp(context_lib->probs[i][n][j]))));
    		}
    		fprintf(fout, "\n");
    	}
    	fputs("//\n", fout);
   	}

    fclose(fout);
}

/* Calculates the log of the probability that profile 'p' emits the counts in
 * a window centered at index 'idx' in 'c'. Note that the normalization factor
 * that is usualy used in multinomial distributions is left out since it
 * cancels out anyway.
 */
int calc_posterior_probabilities_profile(double **profile, int profile_length, int index, CONTEXT_LIB *cs4000, double *posterior_probabilities)
{
	double sum;
	double max = -DBL_MAX;        // Needed to find the maximum number for the log-sum-exp trick
	double posterior_probability; // The posterior probability around a single sequence position

	/* Needed to calculate the posterior probabilites arouns a single sequence position */
	int begin = MAX(0, (index - cs4000->center));
	int end = MIN(profile_length, index + cs4000->center + 1);

	/* Calculate the posterior probabilities */
	for (int i = 0; i < cs4000->size; i++) {
		/* Calculate the posterior probability that the context profile emitted/produced the
		 * observed sequence. The center has more weight than the edges. So it is a simple
		 * calculation of multiplying the probability of the amino acid in a window around index 
		 * being observed in the context profile and multiplying them by the weight of the window
		 */
		posterior_probability = 0.0;
		
		for(int n = begin, k = begin - index + cs4000->center; n < end; ++n, ++k) {
			sum = 0.0;
			for (int a = 0; a < cs4000->nalph; ++a) {
				/* Currenlty we do not substract the backround frequency of the context profile */
				sum += profile[n][a] * cs4000->probs[i][k][a];
			}
		    posterior_probability += cs4000->weights[k] * sum;
		}

		/* Multiply the final probability by the probability of finding the context profile (the prior) */
		posterior_probabilities[i] = cs4000->prior[i] + posterior_probability;

		if (posterior_probabilities[i] > max) {
			max = posterior_probabilities[i];  // needed for log-sum-exp trick
      	}
	}

	// Log-Sum-Exp trick begins here
	sum = 0.0;
	for (int i = 0; i < cs4000->size; i++) {
		sum += exp(posterior_probabilities[i] - max);
	}

	double tmp = max + log(sum);
	for (int i = 0; i < cs4000->size; i++) {
		/* We store the raw probability! */
		posterior_probabilities[i] = exp(posterior_probabilities[i] - tmp);
	}

	return 0;
}

/* Calculates the log of the probability that profile 'p' emits the sequence
 * window centered at index 'idx' in 'seq'. Note that the normalization factor
 * that is usualy used in multinomial distributions is left out since it
 * cancels out anyway.
 */
int calc_posterior_probabilities_sequence(char *sequence, int index, CONTEXT_LIB *cs4000, double *posterior_probabilities)
{
	double max = -DBL_MAX;        // Needed to find the maximum number for the log-sum-exp trick
	double posterior_probability; // The posterior probability around a single sequence position

	/* Needed to calculate the posterior probabilites arouns a single sequence position */
	int begin = MAX(0, (index - cs4000->center));
	int end = MIN(strlen(sequence), index + cs4000->center + 1);

	/* Calculate the posterior probabilities */
	for (int i = 0; i < cs4000->size; i++) {
		/* Calculate the posterior probability that the context profile emitted/produced the
		 * observed sequence. The center has more weight than the edges. So it is a simple
		 * calculation of multiplying the probability of the amino acid in a window around index 
		 * being observed in the context profile and multiplying them by the weight of the window
		 */
		posterior_probability = 0.0;
		for(int n = begin, k = begin - index + cs4000->center; n < end; ++n, ++k) {
			/* Currenlty we do not substract the backround frequency of the context profile */
		    posterior_probability += cs4000->weights[k] * cs4000->probs[i][k][prc_char2int(sequence[n], ALPHABET_PRC)];
		}

		/* Multiply the final probability by the probability of finding the context profile (the prior) */
		posterior_probabilities[i] = cs4000->prior[i] + posterior_probability;

		if (posterior_probabilities[i] > max) {
			max = posterior_probabilities[i];  // needed for log-sum-exp trick
      	}
	}

	// Log-Sum-Exp trick begins here
	double sum = 0.0;
	for (int i = 0; i < cs4000->size; i++) {
		sum += exp(posterior_probabilities[i] - max);
	}

	double tmp = max + log(sum);
	for (int i = 0; i < cs4000->size; i++) {
		/* We store the raw probability! */
		posterior_probabilities[i] = exp(posterior_probabilities[i] - tmp);
	}

	return 0;
}

int init_cs_weights(CONTEXT_LIB *cs, double weight_center, double weight_decay)
{
	if (IS_EVEN(cs->wlen)) {
		fprintf(stderr, "The window length indicated int the context library should be odd, but is %i\n", cs->wlen);
		return 1;
	}

	/* Set the center weight */
	cs->weights[cs->center] = weight_center;

	/* Set the weights using an exponential decay */
	for (int i = 1; i <= cs->center; i++) {
	    double weight = weight_center * pow(weight_decay, i);
	    cs->weights[cs->center - i] = weight;
	    cs->weights[cs->center + i] = weight;
	}

	return 0;
}

double **sequence_to_profile(char *sequence, CONTEXT_LIB *cs4000)
{
	int i, a, k;
	int rank;
	int shape[2];

	/* Allocate space for the final profile */
	shape[0] = strlen(sequence);
	shape[1] = cs4000->nalph;
	rank = 2;
	double **profile = prc_allocarray(rank, shape, sizeof(double));

	/* Allocate space for the pseudocounts */
	shape[0] = cs4000->nalph;
	rank = 1;
	double *pseudocounts = prc_allocarray(rank, shape, sizeof(double));

	/* Allocate space for all the posterior probabilities of single positions in 
	 * sequence for all the context probabilities (size of the library)
	 */
	shape[0] = strlen(sequence);
	shape[1] = cs4000->size;
	rank = 2;
	double **posterior_probabilities = NULL;
	posterior_probabilities = prc_allocarray(rank, shape, sizeof(double));

	/* Initialize weight vector for the calculation of the pseudocounts */
	if (init_cs_weights(cs4000, cs4000->weight_center, cs4000->weight_decay)) {
		fprintf(stderr, "Error initializing weight vector\n");
		return NULL;
	}
	
	/* pointer to a single row containing all the posterior probabilities of all
	 * the context profiles generating a single position of the sequence
	 */
	double *row_of_posterior_probabilities;

	for (i = 0; i < strlen(sequence); i++) {
		row_of_posterior_probabilities = &posterior_probabilities[i][0];

		calc_posterior_probabilities_sequence(sequence, i, cs4000, row_of_posterior_probabilities);
		
		for (k = 0; k < cs4000->size; ++k) {
			for (a = 0; a < cs4000->nalph; ++a) {
				pseudocounts[a] = 0.0;
			}
		}

		for (k = 0; k < cs4000->size; ++k) {
			for (a = 0; a < cs4000->nalph; ++a) {
				/* The amino acid probabilities are weighted against the posterior probabilities */
				/* NOTE: Some posterior probabilities are very small, do we need all of them? */
				pseudocounts[a] += row_of_posterior_probabilities[k] * exp(cs4000->probs[k][cs4000->center][a]);
			}
        }
        /* Finally add the pseudocounts to sequence profile */
        for (a = 0; a < cs4000->nalph; ++a) {
			profile[i][a] = (1.0 - cs4000->pc_admix) * (sequence[i] == ALPHABET_PRC[a] ? 1.0 : 0.0) + cs4000->pc_admix * pseudocounts[a];
		}
	}

	/* We no longer need these */
	free(posterior_probabilities);
	free(pseudocounts);

	/* Return the final profile containing the pseudocounts */
	return profile;
}

unsigned char *profile_to_discrete_alphabet(double **profile, int profile_length, CONTEXT_LIB *cs219)
{
	int i;
	unsigned char *discrete_alphabet = calloc(profile_length, sizeof(unsigned char));

	if (cs219->size > 255) {
		fprintf(stderr, "Discrete alphabet of the CONTEXT_LIB is too big. Only 255 characters allowed\n");
		return NULL;
	}

	/* Initialize the weights for hard clustering of context profiles e.a. translating
	 * a profile to a discrete alphabet using the cs219 context profiles 
	 */
	if (init_cs_weights(cs219, cs219->weight_as, 1.0)) {
		fprintf(stderr, "Could not initialize weights\n");
		return NULL;
	}

	double **posterior_probabilities;
	int rank;
	int shape[2];

	rank = 2;
	shape[0] = profile_length + 1;
	shape[1] = cs219->size;

	posterior_probabilities = prc_allocarray(rank, shape, sizeof(double));
	for (i = 0; i < profile_length; i++) {
		calc_posterior_probabilities_profile(profile, profile_length, i, cs219, &posterior_probabilities[i][0]);
	}

	int index;
	double max;
	for (i = 0; i < profile_length; i++) {
		max = 0.0;
		index = 0;
		for (int n = 0; n < cs219->size; n++) {
			if (posterior_probabilities[i][n] > max) {
				max = posterior_probabilities[i][n];
				index = n;
				if (max > 0.5) {
					break;
				}
			}
		}
		discrete_alphabet[i] = index;
	}
	discrete_alphabet[i] = '\0';

	free(posterior_probabilities);

	return discrete_alphabet;
}

/* Return a new alignment only containing the match columns */
char **extract_match_columns(char **alignment, int nr_sequences, int *match_column)
{
	int i, n, j;
	int rank = 2;
	int shape[2];
	int alilength = strlen(alignment[0]);
	int sum = 0; // The number of match columns
	for (i = 0; i < alilength; i++) {
		sum += match_column[i];
	}

	shape[0] = nr_sequences;
	shape[1] = (int) sum + 1; // the number of match columns + 1
	char **new_alignment;
	new_alignment = prc_allocarray(rank, shape, sizeof(char));

	
	for (i = 0; i < nr_sequences; i++) {
		j = 0;
		for (n = 0; n < strlen(alignment[0]); n++) {
			// If the column is a match column, copy it
			if (match_column[n]) {
				new_alignment[i][j++] = alignment[i][n];
			}
		}
		new_alignment[i][j] = '\0';
	}

	return new_alignment;
}

/* Return a new alignment only containing the insert columns */
char **extract_insert_columns(char **alignment, int nr_sequences, int *match_column)
{
	int i, n, j;
	int rank = 2;
	int shape[2];
	int alilength = strlen(alignment[0]);
	int sum = 0; // The number of match columns
	for (i = 0; i < alilength; i++) {
		sum += match_column[i];
	}
	/* Substract match columns from total colums = insert columns */
	sum = alilength - sum;

	shape[0] = nr_sequences;
	shape[1] = (int) sum + 1; // the number of insert columns + 1
	char **new_alignment;
	new_alignment = prc_allocarray(rank, shape, sizeof(char));

	
	for (i = 0; i < nr_sequences; i++) {
		j = 0;
		for (n = 0; n < strlen(alignment[0]); n++) {
			// If the column is not a match column, copy it
			if (!match_column[n]) {
				new_alignment[i][j++] = alignment[i][n];
			}
		}
	}

	return new_alignment;
}

/* Return a new alignment only containing the deletion columns */
char **extract_deletion_columns(char **alignment, int nr_sequences, int *match_column)
{
	int i, n;
	int rank = 2;
	int shape[2];
	int alilength = strlen(alignment[0]);
	int sum = 0; // The number of match columns
	/* We allocate too much, but it will save us having to loop over the alignment twice */
	for (i = 0; i < alilength; i++) {
		sum += match_column[i];
	}

	shape[0] = nr_sequences;
	shape[1] = (int) sum + 1; // the number of insert columns + 1
	char **new_alignment;
	new_alignment = prc_allocarray(rank, shape, sizeof(char));
	
	int have_deletion;
	int j = 0;
	for (i = 0; i < alilength; i++) {
		have_deletion = 0;
		if (match_column[i]) {
			// Check whether one of the sequences has a gap in the match region e.a. a deletion 
			for (n = 0; n < nr_sequences; n++) {
				if (prc_char2int(alignment[n][i], ANY_OR_GAP) != -1) {
					have_deletion = 1;
				}
			}

			if (have_deletion) {
				for (n = 0; n < nr_sequences; n++) {
					new_alignment[n][j] = alignment[n][i];
				}
				j++;
			}
		}
	}


	return new_alignment;
}

double global_weights_and_diversity(char **alignment, int number_of_sequences, double **global_weights_in, bool neff_sum_pairs)
{
    const double kZero = 1E-10;    // for calculation of entropy
    double *global_weights = NULL; // global sequence weights
    double neff = 0.0f;            // diversity of alignment
    int *nres;                     // number of residues in sequence i
    double *entropy;               // to calculate entropy
    int *ndiffaa;                  // different letters in each column
    int **counts = NULL;           // column counts (excl. ANY)
    int i, n;
    int alilength = strlen(alignment[0]); // Raw length of the alignment

    /* Allocate space for the global weights, entropy, counts, etc */
    global_weights = calloc(number_of_sequences, sizeof(double));
    nres = calloc(number_of_sequences, sizeof(int));
    entropy = calloc(20, sizeof(double));
    ndiffaa = calloc(strlen(alignment[0]), sizeof(int));
    
    int shape[2];
    shape[0] = (int) strlen(alignment[0]);
    shape[1] = 20; // Alphabet size
    int rank = 2;
    counts = prc_allocarray(rank, shape, sizeof(int));

    int temp = 0;
    // Count number of residues in each column
    for (i = 0; i < alilength; i++) {
        for (n = 0; n < number_of_sequences; n++) {
            if ((temp = prc_char2int(alignment[n][i], ALPHABET_PRC)) != -1) {
                /* Add the amino acid observed in the sequence to the column counts */
                counts[i][temp]++;
                /* Increase the amount of residues observed of this sequence */
                nres[n]++;
            }
        }
    }

    // Count number of different residues in each column
    for (i = 0; i < alilength; i++) {
        for (n = 0; n < 20; n++) {
            if (counts[i][n]) {
            	++ndiffaa[i];
            }
        }
        // If there are only gaps and unknown amino acids add a pseudocount of 1
        if (ndiffaa[i] == 0) {
        	ndiffaa[i] = 1;  // col consists of only gaps and ANYs
        }
    }

    // Calculate weights
    for (i = 0; i < alilength; i++) {
        for (n = 0; n < number_of_sequences; n++) {
        	/* If we observed different amino acids and the amino acid is not a gap or unknown
        	 * character, we calculate the global weight of that sequence:
        	 * 1 / number of different amino acids observed at that column * the count at that position * the total sequence length
        	 */
            if (ndiffaa[i] > 0 && ((temp = prc_char2int(alignment[n][i], ALPHABET_PRC)) != -1)) {
                global_weights[n] += 1.0 / (ndiffaa[i] * counts[i][temp] * nres[n]);
            }
        }
    }

    // Normalize
    double sum = prc_sum(global_weights, number_of_sequences);
    for (i = 0; i < number_of_sequences; i++) {
    	global_weights[i] = global_weights[i] / sum;
    }

    // Calculate number of effective sequences
    if (!neff_sum_pairs) {
        for (i = 0; i < alilength; i++) {
        	// Reset entropy
        	for (n = 0; n < 20; n++) {
        		entropy[n] = 0;
        	}

        	/* Calculate the entropy */
            for (n = 0; n < number_of_sequences; n++) {
            	if ((temp = prc_char2int(alignment[n][i], ALPHABET_PRC)) != -1) {
            		entropy[temp] += global_weights[n];
            	}
            }
                
        	// Normalize entropy
        	sum = prc_sum(entropy, 20);
        	for (n = 0; n < 20; n++) {
        		entropy[n] = entropy[n] / sum;
        	}

            for (n = 0; n < 20; n++) {
                if (entropy[n] > kZero) {
                	neff -= entropy[n] * log2(entropy[n]);
                }
            }
        }
        neff = pow(2.0, neff / alilength);
    } else {
        for (i = 0; i < alilength; i++) {
        	// Reset entropy
        	for (n = 0; n < 20; n++) {
        		entropy[n] = 0;
        	}

            for (n = 0; n < number_of_sequences; n++) {
                if ((temp = prc_char2int(alignment[n][i], ALPHABET_PRC)) != -1) {
                	entropy[temp] += global_weights[n];
                }
            }

        	// Normalize entropy
        	sum = prc_sum(entropy, 20);
        	for (n = 0; n < 20; n++) {
        		entropy[n] = entropy[n] / sum;
        	}

            double ni = number_of_sequences + 1.0;

            for (n = 0; n < 20; n++) {
                ni -= prc_sqr(entropy[n]) * number_of_sequences;
            }

            neff += ni;
        }
        neff /= alilength;
    }

    free(nres);
    free(entropy);
    free(ndiffaa);
    free(counts);

    *global_weights_in = global_weights;

    return neff;
}

// Calculates position-dependent sequence weights (pswm, alignment length rows, number of sequence columns) and 
// returns the number of effective sequences (alignment length long)
double *position_specific_weights_and_diversity(char **alignment, int number_of_sequences, double ***pswm_in)
{
//     // Maximal fraction of sequences with an endgap (Currently we do not deal with ENDGAPs!)
//     const double kMaxEndgapFraction = 0.1;
	int i, n, j, k = 0;
	int rank = 2;
	int shape[2];

    const double kZero = 1E-10; // Zero value for calculation of entropy
    int min_cols = 10;    // Minimum number of columns in subalignments
    int alilength = strlen(alignment[0]); // Raw length of the alignment

    // Return values
    double *neff;  // Diversity of subalignment i
    double **pswm; // The position dependent sequence weights of seq n in column i

    // Helper variables
    int ncoli = 0;        // number of columns j that contribute to neff[i]
    int nseqi = 0;        // number of sequences in subalignment i
    int ndiff = 0;        // number of different alphabet letters
    double sum = 0;       // used for normalization
    double *entropy;               // to calculate entropy
    bool change = false;  // has the set of sequences in subalignment changed?

    // Number of seqs with some residue in column i AND a at position j
    int **nres;
    shape[0] = (int) strlen(alignment[0]);
    shape[1] = strlen(ALPHABET_PRC); // Alphabet size including gaps and any character
    rank = 2;
    nres = prc_allocarray(rank, shape, sizeof(int));

    // Weight of sequence k in column i, calculated from subalignment i
    double *weight;
    weight = calloc(number_of_sequences, sizeof(double));

    /* Allocate space for the weights, entropy, counts, etc */
    neff = calloc(alilength, sizeof(double));
    entropy = calloc(20, sizeof(double));
    
    shape[0] = (int) strlen(alignment[0]);
    shape[1] = number_of_sequences; // Alphabet size
    rank = 2;
    pswm = prc_allocarray(rank, shape, sizeof(double));

    double *global_weights;
	global_weights_and_diversity(alignment, number_of_sequences, &global_weights, 1);

    int temp; // holds the index of the amino acid
    for (i = 0; i < alilength; i++) {
        change = false;
        for (n = 0; n < number_of_sequences; ++n) {
        	/* increase the number of sequences at start or if the previous column 
        	 * has a gap and the current an amino acid (not a gap or unknown)
        	 */
			if ((i == 0 || (prc_char2int(alignment[n][i - 1], ANY_OR_GAP) > -1)) 
            	&& ((prc_char2int(alignment[n][i], ANY_OR_GAP)) == -1)) {
                change = true;
                ++nseqi;
                for (j = 0; j < alilength; j++) {
                	temp = prc_char2int(alignment[n][j], ALPHABET_PRC);
                	if (temp > -1) {
                    	++nres[j][temp];
                	}
                }
         	/* decrease the number of sequences if the previous column 
        	 * has is not a gap or unknown and the current is a gap
        	 */               
            } else if (i > 0 && (prc_char2int(alignment[n][i - 1], ANY_OR_GAP) == -1)
            	&& (prc_char2int(alignment[n][i], ANY_OR_GAP) != -1)) {    
                change = true;
                --nseqi;
                for (j = 0; j < alilength; ++j) {
                	temp = prc_char2int(alignment[n][j], ALPHABET_PRC);
                	if (temp > -1) {
                    	--nres[j][temp];
                   	}
                }
            }
        }

        if (change) {  // set of sequences in subalignment has changed
            ncoli = 0;
            
            //Reset the weights
            for (n = 0; n < number_of_sequences; n++) {
            	weight[n] = 0;
            }

            for (j = 0; j < alilength; j++) {
            	// Currently we are not handling ENDGAPs!
                //if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
                ndiff = 0;
                for (n = 0; n < 20; n++) {
                	if (nres[j][n]) {
                		++ndiff;
                	}
                }

                if (ndiff == 0) {
                	continue;
                }

                ++ncoli;
                
                for (n = 0; n < number_of_sequences; n++) {
                	// Make sure both alignment positions contain an amino acid
                    if ((prc_char2int(alignment[n][i], ANY_OR_GAP) == -1) && 
                    	(prc_char2int(alignment[n][j], ANY_OR_GAP) == -1)) {
                    	temp = prc_char2int(alignment[n][j], ALPHABET_PRC);
                    	weight[n] += 1.0 / (double)((nres[j][temp] * ndiff));
                    }
                }
            }  // for j over ncols

        	// Normalize weights
        	sum = prc_sum(weight, number_of_sequences);
        	for (n = 0; n < number_of_sequences; n++) {
        		weight[n] = weight[n] / sum;
        	}

            if (ncoli < min_cols) { // number of columns in subalignment insufficient?
                for (k = 0; k < number_of_sequences; k++) {
                    if ((prc_char2int(alignment[k][i], ANY_OR_GAP) == -1)) {
                        weight[k] = global_weights[k];
                    } else {
                        weight[k] = 0.0f;
                    }
                }
            }

            neff[i] = 0.0f;
            for (n = 0; n < alilength; n++) {
                //if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
               
                //Reset entropy
            	for (j = 0; j < 20; j++) {
            		entropy[j] = 0.0;
            	}

                for (k = 0; k < number_of_sequences; k++) {
                    if ((prc_char2int(alignment[k][i], ANY_OR_GAP) == -1) && 
                    	(prc_char2int(alignment[k][n], ANY_OR_GAP) == -1)) {
                    	temp = prc_char2int(alignment[k][n], ALPHABET_PRC);
                        entropy[temp] += weight[k];
                    }
                }

                //Normalize entropy
                sum = prc_sum(entropy, 20);
            	for (j = 0; j < 20; j++) {
            		entropy[j] = entropy[j] / sum;
            	}

                for (j = 0; j < 20; j++) {
                    if (entropy[j] > kZero) {
                    	neff[i] -= entropy[j] * log2(entropy[j]);
                   	}
                }
            }  // for j over ncols

            neff[i] = (ncoli > 0 ? pow(2.0, neff[i] / ncoli) : 1.0f);

        } else {  // set of sequences in subalignment has NOT changed
            neff[i] = (i == 0 ? 0.0 : neff[i-1]);
        }

        for (k = 0; k < number_of_sequences; k++) {
        	pswm[i][k] = weight[k];
        }
    }  // for i over ncols

    // Clean-up
    free(nres);
    free(weight);
    free(entropy);
    free(global_weights);
    
    *pswm_in = pswm;

	return neff;
}

/* This function returns a int array indicating whether the column (location in the sequence) is a match column or
 * not. This is determined using global weights and a gap threshold (either between 0.0 and 1.0 as a fraction of
 * as a number between 1.1 and 100 as a percentage).
 */
int *assign_match_columns_by_gap_rule(char **alignment, int number_of_sequences, double gap_threshold)
{
	double neff;
	double *global_weights;
	int *match_columns = NULL;
	int alilength = strlen(alignment[0]);
	int i, n;

	match_columns = calloc(alilength + 1, sizeof(int));

	neff = global_weights_and_diversity(alignment, number_of_sequences, &global_weights, 1);

	int temp = 0;
	for (i = 0; i < alilength; i++) {
		double gap = 0.0f;
		double residue = 0.0f;

		for (n = 0; n < number_of_sequences; n++) {
			if ((temp = prc_char2int(alignment[n][i], ALPHABET_PRC)) != -1) {
            	residue += global_weights[n];
            } else {
            	// In HHsearch the ENDGAPs are ignored, should we do the same?
            	gap += global_weights[n];
            }
		}

        if (gap_threshold > 1.0) { // interpret as number between 1 and 100
            double percent_gaps = 100.0 * gap / (residue + gap);
            match_columns[i] = (percent_gaps <= gap_threshold);
        } else { // interpret as decimal number
            double frac_res = residue / (residue + gap);
            match_columns[i] = (frac_res > gap_threshold);
        }
	}

	match_columns[i] = 999; // Gard for terminal insertion states, something high that we don't use

	free(global_weights);

	return match_columns;
}

/* Calculate position-dependent sequence weights (pswm) and the number of effective
 * sequences (neff). Otherwise many similar sequences will bias the informational content of the 
 * HMM too much. This is done according to Henikoff 1999. This is then used to adjust the amino acid count
 */ 
double **alignment_to_count_profile(char **alignment, int number_of_sequences, double **neff_in, bool pos_weight, bool neff_sum_pairs)
{
	int i, n;
	int alilength = strlen(alignment[0]);
	int rank = 2;
	int shape[2];
	shape[0] = alilength;
	shape[1] = 20; // The alphabet size 
	double **counts = NULL;
	counts = prc_allocarray(rank, shape, sizeof(double));
	int temp;
	// Add counts and neff from alignment to count profile
	if (pos_weight) {
		double **pswm; // The position specific weight matrix
		double *neff;  // The number of effective sequences

		neff = position_specific_weights_and_diversity(alignment, number_of_sequences, &pswm);

		/* Return neff, needed for pseudocounts */
		*neff_in = neff;
		
		for (i = 0; i < alilength; i++) {
			for (n = 0; n < number_of_sequences; n++) {
				if (prc_char2int(alignment[n][i], ANY_OR_GAP) == -1) {
					temp = prc_char2int(alignment[n][i], ALPHABET_PRC);
					counts[i][temp] += pswm[i][n] * neff[i];
				}
			}
		}
	} else {  // use faster global sequence weights
		double *global_weights;
		double global_neff;
		global_neff = global_weights_and_diversity(alignment, number_of_sequences, &global_weights, neff_sum_pairs);
		for (i = 0; i < alilength; i++) {
			for (n = 0; n < number_of_sequences; n++) {
				if (prc_char2int(alignment[n][i], ANY_OR_GAP) == -1) {
					temp = prc_char2int(alignment[n][i], ALPHABET_PRC);
					counts[i][temp] += global_weights[n] * global_neff;
				}
			}
		}		
	}

	return counts;
}

// AddToProfile(const CountProfile<Abc>& cp,
//                                             const Admix& pca,
//                                             Profile<Abc>& p) const {
//     assert_eq(cp.counts.length(), p.length());
//     LOG(INFO) << "Adding library pseudocounts to profile ...";

//     Matrix<double> pp(cp.counts.length(), lib_.size(), 0.0);  // posterior probs
//     Vector<double> pc(Abc::kSize, 0.0);                       // pseudocount vector P(a|X_i)

//     // Calculate and add pseudocounts for each sequence window X_i separately
//     for (size_t i = 0; i < cp.counts.length(); ++i) {
//         double* ppi = &pp[i][0];
//         // Calculate posterior probability of state k given sequence window around 'i'
//         CalculatePosteriorProbs(lib_, emission_, cp, i, ppi);
//         // Calculate pseudocount vector P(a|X_i)
//         Assign(pc, 0.0);
//         for (size_t k = 0; k < lib_.size(); ++k)
//             for(size_t a = 0; a < Abc::kSize; ++a)
//                 pc[a] += ppi[k] * lib_[k].pc[a];
//         // FIXME: is normalization here really needed?
//         Normalize(&pc[0], Abc::kSize);
//         // Add pseudocounts to profile
//         double tau = pca(cp.neff[i]);
//         for(size_t a = 0; a < Abc::kSize; ++a)
//             p[i][a] = (1.0 - tau) * cp.counts[i][a] / cp.neff[i] + tau * pc[a];
//     }
// }

double **count_profile_to_profile(double **count_profile, int profile_length, double *neff, CONTEXT_LIB *cs4000)
{
	int i, k, a;
	int rank = 2;
	int shape[2];
	shape[0] = profile_length;
	shape[1] = 20; // alphabet size
	double **profile = NULL;
	profile = prc_allocarray(rank, shape, sizeof(double));

	/* Allocate space for the pseudocounts */
	shape[0] = cs4000->nalph;
	rank = 1;
	double *pseudocounts = prc_allocarray(rank, shape, sizeof(double));

	/* Initialize the weights for soft clustering of context profiles e.a. 
	 * adding pseudocounts and NOT translating to disrete alphabet
	 */
	if (init_cs_weights(cs4000, cs4000->weight_center, cs4000->weight_decay)) {
		fprintf(stderr, "Could not initialize weights\n");
		return NULL;
	}

	double **posterior_probabilities;
	rank = 2;
	shape[0] = profile_length + 1;
	shape[1] = cs4000->size;
	
	/* Allocate space for the posterior probabilities */
	posterior_probabilities = prc_allocarray(rank, shape, sizeof(double));

	/* pointer to a single row containing all the posterior probabilities of all
	 * the context profiles generating a single position of the sequence
	 */
	double *row_of_posterior_probabilities;
	for (i = 0; i < profile_length; i++) {
		row_of_posterior_probabilities = &posterior_probabilities[i][0];
		calc_posterior_probabilities_profile(count_profile, profile_length, i, cs4000, row_of_posterior_probabilities);
		
		for (k = 0; k < cs4000->size; ++k) {
			for (a = 0; a < cs4000->nalph; ++a) {
				pseudocounts[a] = 0.0;
			}
		}

		for (k = 0; k < cs4000->size; ++k) {
			for (a = 0; a < cs4000->nalph; ++a) {
				/* The amino acid probabilities are weighted against the posterior probabilities */
				/* NOTE: Some posterior probabilities are very small, do we need all of them? */
				pseudocounts[a] += row_of_posterior_probabilities[k] * exp(cs4000->probs[k][cs4000->center][a]);
			}
        }
        /* Finally add the pseudocounts to sequence profile */
        double pca = cs4000->pc_admix; // = 0.3
        double pcb = cs4000->pc_ali; // = 4.0
		double tau = MIN(1.0, pca * (pcb + 1.0) / (pcb + neff[i]));
		//p[i][a] = (1.0 - tau) * cp.counts[i][a] / cp.neff[i] + tau * pc[a];
        for (a = 0; a < cs4000->nalph; ++a) {
			profile[i][a] = (1.0 - tau) * count_profile[i][a] / neff[i] + tau * pseudocounts[a];
		}

	}

	free(posterior_probabilities);
	free(pseudocounts);

	return profile;
}

int make_profile_transitions_and_neff(char **alignment, int nr_sequences, double ***profile_in, double ***transitions_in, double ***NeffMID_in)
{
    int i, n;
    int *match_column = NULL;
    double gap_threshold = 0.5;
    /* Determine which locations in the alignment are gaps using 
     * global entropy weights according to Henikoff&Henikoff '94 
     */
    match_column = assign_match_columns_by_gap_rule(alignment, nr_sequences, gap_threshold);
    
    /* The new HMM will be as big as the match columns */
    int hmmlength = 0;
    for (i = 0; i < strlen(alignment[0]); i++) {
    	hmmlength += match_column[i];
    }

    /* Now that we know what positions are gaps, make a new alignment without them. Only
     * keep the match positions, so we can count them and add a pseudocount later.
     */
    char **compressed_alignment = extract_match_columns(alignment, nr_sequences, match_column);

    /* Use the newly formed alignment to detemine the position specific weight matrix (or global weights) and
     * use that to return a count profile (the amount of aa per position) corrected for the sequence weight
     */
    double **count_profile = NULL;
    bool pos_weight = true;
    bool neff_sum_pairs = true;
    double *neff;
    count_profile = alignment_to_count_profile(compressed_alignment, nr_sequences, &neff, pos_weight, neff_sum_pairs);

    /* Now we can finally continue with the next step and add the context specific pseudocount.
     * We will start by loading the context library.
     */
 	CONTEXT_LIB *cs4000;
	cs4000 = read_context_library("cs4000.lib");
	//write_context_library("tempcontext_dataprc.lib", cs4000);

	// add context sensitive pseudocounts
	double **profile = NULL; // contains the probabilities of the amino acids in a sequence length (rows) by number of amino acids matrix (columns);
	int profile_length = strlen(compressed_alignment[0]);
	profile = count_profile_to_profile(count_profile, profile_length, neff, cs4000);
	*profile_in = profile;

	/* Transitions and effective sequences are calculated using position specific weight matrix 
	 * for match states, but global weights for deletions and insertions 
	 */
	// Get the global weights
	double *global_weights;
	double global_neff;
	global_neff = global_weights_and_diversity(compressed_alignment, nr_sequences, &global_weights, 1);

	const int INSERT_STATE = 0;
	const int MATCH_STATE = 1;
	const int DELETION_STATE = 2;
	const int NUM_STATES = 3;

	// Make the NeffMID matrix
	int rank = 2;
	int shape[2];
	shape[0] = NUM_STATES;
	shape[1] = strlen(alignment[0]);
	double **NeffMID;
	NeffMID = prc_allocarray(rank, shape, sizeof(double));

	*NeffMID_in = NeffMID;

	const int transM2M = 0;
	const int transM2I = 1;
	const int transM2D = 2;
	const int transI2M = 3;
	const int transI2I = 4;
	const int transD2M = 5;
	const int transD2D = 6;
	const int TRANSITIONS = 7;

	shape[0] = strlen(alignment[0]);
	shape[1] = TRANSITIONS;
	double **transitions;
	transitions = prc_allocarray(rank, shape, sizeof(double));

	*transitions_in = transitions;

	double **pswm; // The position specific weight matrix
	free(neff);

	neff = position_specific_weights_and_diversity(compressed_alignment, nr_sequences, &pswm);

    float Nlim = fmax(10.0, global_neff + 1.0);    // limiting Neff
    float scale = log2((Nlim - global_neff) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos


	int number_of_insert_states = 0;
	int incols = 0;
	int dncols = 0;
	double iweight = -1.0 / nr_sequences;
	double dweight = -1.0 / nr_sequences;
	double sum = 0.0;
	int *amount_amino_acids_in_insert_state_sequences;
	amount_amino_acids_in_insert_state_sequences = calloc(nr_sequences, sizeof(int));
	for (i = 0; i < strlen(alignment[0]); i++) {
		dncols = 0;
		dweight = -1.0 / nr_sequences;
		sum = 0.0;

		for (n = 0; n < nr_sequences; n++) {
			if (match_column[i] == MATCH_STATE) {
				if (prc_char2int(alignment[n][i], ANY_OR_GAP) == -1) { // current state is M
					/* Deletion to match and match to match can span insertion gaps */
					int found = 0;
					if (match_column[i + 1] == INSERT_STATE) {
						int temp = 1;
						while (match_column[i + temp] == INSERT_STATE) {
							if (prc_char2int(alignment[n][i + temp], GAPS) == -1) {
								if (!found) {
								// We have a match to insert
								transitions[i - number_of_insert_states][transM2I] += pswm[i - number_of_insert_states][n];
								}
							}
							temp++;
						}
						if (!found) {
							if ((prc_char2int(alignment[n][i + temp], GAPS) > -1) && (match_column[i + temp] == MATCH_STATE)) { //next state is D
								transitions[i - number_of_insert_states][transM2D] += pswm[i - number_of_insert_states][n];
							} else if (prc_char2int(alignment[n][i + temp], ALPHABET_PRC) > -1) { //next state is M
								transitions[i  - number_of_insert_states][transM2M] += pswm[i - number_of_insert_states][n];
							}
						}
					} else {
						if ((prc_char2int(alignment[n][i + 1], GAPS) > -1) && (match_column[i + 1] == MATCH_STATE)) { //next state is D
							transitions[i - number_of_insert_states][transM2D] += pswm[i - number_of_insert_states][n];
						} else if (prc_char2int(alignment[n][i + 1], ALPHABET_PRC) > -1) { //next state is M
							transitions[i  - number_of_insert_states][transM2M] += pswm[i - number_of_insert_states][n];
						} else if (prc_char2int(alignment[n][i + 1], GAPS) > -1) { //next state is I
							transitions[i - number_of_insert_states][transM2I] += pswm[i - number_of_insert_states][n];
						}
					}
				} else if (prc_char2int(alignment[n][i], GAPS) > -1) { // current state is D
					dncols++;
					dweight += global_weights[n];
					/* Deletion to match and match to match can span insertion gaps */
					if (match_column[i + 1] == INSERT_STATE) {
						int temp = 1;
						while (match_column[i + temp] == INSERT_STATE) {
							temp++;
						}
						if (prc_char2int(alignment[n][i + temp], GAPS) > -1) { //next state after insertion is D
							transitions[i - number_of_insert_states][transD2D] += global_weights[n];
						} else if (prc_char2int(alignment[n][i + temp], ALPHABET_PRC) > -1) { //next state after insertion is M
							transitions[i - number_of_insert_states][transD2M] += global_weights[n];
						}
					} else {
						if (prc_char2int(alignment[n][i + 1], GAPS) > -1) { //next state is D
							transitions[i - number_of_insert_states][transD2D] += global_weights[n];
						} else if (prc_char2int(alignment[n][i + 1], ALPHABET_PRC) > -1) { //next state is M
							transitions[i - number_of_insert_states][transD2M] += global_weights[n];
						}
					}
				}
			} else if (match_column[i] == INSERT_STATE) {
				if (prc_char2int(alignment[n][i], GAPS) == -1) { // we only use those sequences with and amino acid to determine the entropy
					amount_amino_acids_in_insert_state_sequences[n]++;
				}

				if (prc_char2int(alignment[n][i], GAPS) > -1) {  //cuurent state is I
					if ((prc_char2int(alignment[n][i + 1], GAPS) == -1)  && (match_column[i + 1] == MATCH_STATE)) {  //next state is M
						transitions[i - number_of_insert_states][transI2M] += global_weights[n];
					} else if (prc_char2int(alignment[n][i + 1], GAPS) > -1) { //next state is I {
						transitions[i - number_of_insert_states][transI2I] += global_weights[n];
					}					
				}
			}
		}

		// /* Number of effective sequences */
		if (dncols > 0) {
			if (dweight < 0.0) {
				NeffMID[DELETION_STATE][i - number_of_insert_states] = 1.0;
			} else {
            	NeffMID[DELETION_STATE][i - number_of_insert_states] = Nlim - (Nlim-1.0)*pow(2.0, scale*dweight);
           	}
		}

		// Leaving the insert state
		if ((match_column[i + 1] == MATCH_STATE) && (match_column[i + 1] == MATCH_STATE)) {
			for (n = 0; n < nr_sequences; n++) {
				if (amount_amino_acids_in_insert_state_sequences[n]) {
					incols++;
					iweight += global_weights[n];
				}

				// Reset counts
				amount_amino_acids_in_insert_state_sequences[n] = 0;
			}

			if (incols > 0) {
				if (iweight < 0.0) {
					NeffMID[INSERT_STATE][i - number_of_insert_states] = 1.0;
				} else {
					NeffMID[INSERT_STATE][i - number_of_insert_states] = Nlim - (Nlim-1.0)*pow(2.0, scale*iweight);
				}
			}
			incols = 0;
			iweight = -1.0 / nr_sequences;
		}


		NeffMID[MATCH_STATE][i - number_of_insert_states] = neff[i - number_of_insert_states];

		// for all sequences
		if (match_column[i] == INSERT_STATE) {
			number_of_insert_states++;
		}
	}

	for (i = 0; i < strlen(alignment[0]); i++) {
		/* Normalize transitions */
		// Match states
		sum = transitions[i][transM2D] + transitions[i][transM2I] + transitions[i][transM2M];
		if (sum) {
			transitions[i][transM2D] /= sum;
			transitions[i][transM2I] /= sum;
			transitions[i][transM2M] /= sum;
		}
		// Insert states
		sum = transitions[i][transI2M] + transitions[i][transI2I];
		if (sum) {
			transitions[i][transI2M] /= sum;
			transitions[i][transI2I] /= sum;
		}
		// Deletion states
		sum = transitions[i][transD2M] + transitions[i][transD2D];
		if (sum) {
			transitions[i][transD2M] /= sum;
			transitions[i][transD2D] /= sum;
		}

// TODO Double check this!!!
		// printf("---%i---\n", i);
		// printf("M2D: %lf\n", transitions[i][transM2D]);
		// printf("M2I: %lf\n", transitions[i][transM2I]);
		// printf("M2M: %lf\n", transitions[i][transM2M]);
		// printf("I2M: %lf\n", transitions[i][transI2M]);
		// printf("I2I: %lf\n", transitions[i][transI2I]);
		// printf("D2M: %lf\n", transitions[i][transD2M]);
		// printf("D2D: %lf\n", transitions[i][transD2D]);

		// printf("NeffM: %lf\n",NeffMID[MATCH_STATE][i]);
		// printf("NeffI: %lf\n",NeffMID[INSERT_STATE][i]);
		// printf("NeffD: %lf\n",NeffMID[DELETION_STATE][i]);
	}
	*NeffMID_in = NeffMID;

	free(amount_amino_acids_in_insert_state_sequences);
	free(match_column);
	free(compressed_alignment);
	free(neff);
	free(count_profile);
	free_context_library(cs4000);

    return hmmlength;

}
