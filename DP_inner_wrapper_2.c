/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

    DP_inner_wrapper_2.c: generates different code based on p->MM_function

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2004-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

if (dp->p->MM_function == DOT1)
{
#define MM_FUNCTION DOT1
#include "DP_inner_wrapper_1.c"
#undef  MM_FUNCTION
} else if (dp->p->MM_function == DOT2)
{
#define MM_FUNCTION DOT2
#include "DP_inner_wrapper_1.c"
#undef  MM_FUNCTION
}
