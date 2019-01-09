/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                   PRCx, the profile comparer eXtended, version 1.0.0

      DP_inner_wrapper_1.c: hand-coded extension of auto_DP_inner_loop.c

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Copyright (C) 2004-5 Martin Madera and MRC LMB, Cambridge, UK
   Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
   All Rights Reserved

   This source code is distributed under the terms of the GNU General Public
   License. See the files COPYING and LICENSE for details.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if   MM_FUNCTION == DOT1
#define M1 hmm1->Pmat[m1]
#define M2 hmm2->Pmat[m2]
#define I1 hmm1->Pins[m1]
#define I2 hmm2->Pins[m2]
#define SCORE(x,y)  score_dot1(x, y, null_score)
#elif MM_FUNCTION == DOT2
#define M1 hmm1->Smat[m1]
#define M2 hmm2->Smat[m2]
#define I1 hmm1->Sins[m1]
#define I2 hmm2->Sins[m2]
#define SCORE(x,y)  score_dot2(x, y)
#endif

for (m1 = m1_11_row; m1_end_cond; m1_diff)
{
	for (m2 = m2_11_col; m2_end_cond; m2_diff) {
#if ALGORITHM != OPT_ACCUR
#if X_STATUS == X_NOT_SET

		X[m1][m2][xMM] = SCORE(M1, M2);
		X[m1][m2][xMI] = SCORE(M1, I2);
		X[m1][m2][xIM] = SCORE(I1, M2);
#if N_X_STATES == 4
		X[m1][m2][xII] = SCORE(I1, I2);
#endif

#elif X_STATUS == X_LIN2LOG

		X[m1][m2][xMM] = LOG(X[m1][m2][xMM]);
		X[m1][m2][xMI] = LOG(X[m1][m2][xMI]);
		X[m1][m2][xIM] = LOG(X[m1][m2][xIM]);
#if N_X_STATES == 4
		X[m1][m2][xII] = LOG(X[m1][m2][xII]);
#endif

#endif // #if X_STATUS == ...
#endif // #if ALGORITHM != OPT_ACCUR

		StBBMM = St1[m1][tBM] x St2[m2][tBM];
		StMMEE = St1[m1][tME] x St2[m2][tME];

#include "auto_DP_inner_loop.c"
	};
};

#undef M1
#undef M2
#undef I1
#undef I2
#undef SCORE

#undef xm1_M
#undef xm1_I
#undef xm2_M
#undef xm2_I
