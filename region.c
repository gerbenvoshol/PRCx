/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

	    region.c: two functions for manipulating struct REGION

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <string.h>
#include "prc.h"


// N.B. the order is the same as for cp: src -> dst
//
void copy_region(REGION *src, REGION *dst)
{
	memcpy((void*)dst, (void*)src, sizeof(REGION));
}

void set_region(REGION *r, int start1, int end1, int start2, int end2)
{
	r->start1 = start1;
	r->end1   = end1;
	r->start2 = start2;
	r->end2   = end2;
}
