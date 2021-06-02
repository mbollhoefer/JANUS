/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>


void PRINTMATRIX(SPARSEMATRIX *A) {

  integer i,j,n=A->nc;

  for (i=0; i<n; i++) {
      printf("column %3ld\n",i);
      for (j=0; j<A->ncol[i]; j++) {
	printf("%12ld",A->rowind[i][j]);
      }
      printf("\n");
      if (A->val!=NULL) {
	 if (A->val[i]!=NULL) {
	    for (j=0; j<A->ncol[i]; j++) {
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_
	        printf("%12.4le",A->val[i][j].r);
#else
		printf("%12.4le",A->val[i][j]);
#endif
	    }
	    printf("\n");
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_
	    for (j=0; j<A->ncol[i]; j++) {
	        printf("%12.4le",A->val[i][j].i);
	    }
	    printf("\n");
#endif
	    fflush(stdout);
	 }
      }
  }
} // end PRINTMATRIX
