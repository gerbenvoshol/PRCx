/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

	      matches.c: functions for manipulating struct MATCH

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2004-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <stdlib.h>
#include "prc.h"


// free a given match and return its next
//
MATCH* free_match(MATCH *match)
{
	MATCH *next = match->next;

	if (match->free_hmm2i) {
		free_HMMinfo(match->hmm2i);
	}

	free_2D_array(match->trace);
	free_unless_null(match->alignment);
	free_unless_null(match);

	return next;
}

// free all members of a list of matches
//
void free_matches(MATCH *match)
{
	while (match) {
		match = free_match(match);
	}
}

// count the number of matches in a list
//
int count_matches(MATCH *list)
{
	MATCH *match;
	int   n = 0;

	for (match = list; match; match = match->next) {
		n++;
	}

	return n;
}

// join two lists & return the joint list
//
// the shorter list should be passed first
//
MATCH* join_lists(MATCH *list1, MATCH *list2)
{
	MATCH *match;

	if (list1 == NULL) {
		return list2;
	}

	match = list1;
	while (match->next) {
		match = match->next;
	}

	match->next = list2;

	return list1;
}

// compare two matches based on E-value or reverse score
//
int comp_matches_E_value(const void *match1, const void *match2)
{
	double E1 = (*((MATCH**)match1))->E_value;
	double E2 = (*((MATCH**)match2))->E_value;

	if (E1 < E2) {
		return -1;
	} else if (E1 == E2) {
		return 0;
	}

	return 1;
}
int comp_matches_reverse(const void *match1, const void *match2)
{
	double r1 = (*((MATCH**)match1))->reverse;
	double r2 = (*((MATCH**)match2))->reverse;

	if (r1 < r2) {
		return 1;
	} else if (r1 == r2) {
		return 0;
	}

	return -1;
}

// take a list of matches, copy them to p->matches, and sort the array
//
// returns 0 if all OK, 1 otherwise
//
int copy_sort_matches(PARAMS *p, MATCH *list)
{
	MATCH *match;
	long  int i = 0;

	p->n_matches = count_matches(list);

	if (p->n_matches == 0) {
		return 1;
	}

	malloc_1D_array(p->matches, MATCH*, p->n_matches);

	for (match = list; match; match = match->next) {
		p->matches[i++] = match;
	}

	if (p->E_values)
		qsort((void*) p->matches, p->n_matches, sizeof(MATCH*),
		      comp_matches_E_value);
	else
		qsort((void*) p->matches, p->n_matches, sizeof(MATCH*),
		      comp_matches_reverse);

	return 0;
}
