/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   PRCx, the profile comparer eXtended, version 1.0.0

	 common_main.c: routines common to prc.c and convert_to_prc.c

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2004-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include "prc.h"

char version[100];

// print out the header
//
void print_header(FILE *stream, char *init)
{
	fprintf(stream,

	        "%s%s\n"
	        "%s" COPYRIGHT "\n"
	        "%s" COPYRIGHT2 "\n"
	        "%s" LICENSE   "\n"
	        "%s\n",

	        init, version,
	        init,
	        init,
	        init,
	        init);
};

// print out the help message & die
//
void print_help_die()
{
	printf("%s", usage);
	exit(1);
}

// print out an erorr message and die
//
void arg_error(char *fmt, ...)
{
	va_list list;

	va_start(list, fmt);

	fprintf(stderr, "\n");
	vfprintf(stderr, fmt, list);
	fprintf(stderr, "\n\n\n(To see the help message, "
	        "run the program with no arguments.)\n\n");

	va_end(list);

	exit(1);
}
