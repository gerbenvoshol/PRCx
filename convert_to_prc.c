/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   PRCx, the profile comparer eXtended, version 1.0.0

	convert_to_prc.c: wrapper for converter to the PRC file format

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2003-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prc.h"


static char usage[] =
    "Usage: convert_to_prc [-human] [-name] [<HMM name>] <input file> [<output file>]\n"
    "\n"
    "Recognized input file extensions:\n"
    "\n"
    KNOWN_EXTENSIONS
    "\n"
    "The default output format is the current PRC binary format. When\n"
    "the '-human' option is used, the output will be a human-readable\n"
    "ASCII format similar to SAM. The '-human' format cannot be read\n"
    "by PRC, but it may be useful for debugging and experimentation.\n"
    "The name of the HMM can optionally be given usig the -name option.\n"
    "\n"
    "When no output file is specified, the output is to STDOUT.\n"
    "\n";


// need to have usage[] for common.c
#include "common.c"


int main(int argc, char **argv)
{
	FILE  *file;
	int   argi, human = 0, name_index = 0; // Defaults
	HMM   *hmm;

	sprintf(version, "CONVERT-TO-PRC " VERSION " (%s), compiled on %s",
#if   PROF_HMM_TRANS==PLAN9
	        "PLAN9"
#elif PROF_HMM_TRANS==PLAN7
	        "PLAN7"
#endif
	        , __DATE__);

	// print out the header
	print_header(stderr, "");

	// sort out the options
	for (argi = 1; argi < argc; argi++) {
		// onto the filenames ...
		if (argv[argi][0] != '-') {
			break;
		}
		if (strncmp(argv[argi], "-human", 7) == 0) {
			human = 1;
		} else if (strncmp(argv[argi], "-name", 6) == 0) {
			name_index = ++argi;
		} else if ((strncmp(argv[argi], "-h", 3) == 0) ||
		           (strncmp(argv[argi], "--help", 7) == 0)) {
			print_help_die();
		} else {
			arg_error("Error: Unknown option '%s'!", argv[argi]);
		};
	};

	// check that have one or two filenames
	if (argc - argi == 0) {
		print_help_die();
	} else if (argc - argi > 2) {
		arg_error("Error: Incorrect number of arguments!");
	}


	hmm = read_HMM(argv[argi++]);

	/* If we want to set the name */
	if (name_index) {
		if (hmm->i->name) {
			free(hmm->i->name);
		}
		hmm->i->name = strdup(argv[name_index]);
	}

	if (argc - argi == 1) {
		open_file_or_die(file, argv[argi], "w");
	} else {
		file = stdout;
	};

	if (human) {
		/* If we don't have a discrete alphabet yet, make one */
		if (hmm->discrete_alphabet == NULL) {
			CONTEXT_LIB *cs219 = NULL;
			if ((cs219 = read_context_library("cs219.lib")) == NULL) {
				fprintf(stderr, "Error reading context library from: %s\n", "cs219.lib");
				exit(1);
			}

			hmm->discrete_alphabet = profile_to_discrete_alphabet((double **)&hmm->Pmat[2], hmm->i->M, cs219);
			if (hmm->discrete_alphabet == NULL) {
				fprintf(stderr, "Error converting HMM emissions to discrete alphabet\n");
				exit(1);			
			}
		}
		print_HMM(file, hmm);
	} else {
		// Old file format without discrete alphabet
		//write_HMM_PRC_binary(file, hmm);
		write_HMM_PRC_binary_with_discrete_alphabet(file, hmm);
	}

	if (argc - argi == 1) {
		fclose(file);
	}

	free_HMM(hmm);

	return 0;
}
