/* $Id: symmatvec.c 3608 2017-09-05 07:33:23Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)
#define TRANSA        "t"

#ifdef _DOUBLE_REAL_
#define MYSYMMATVEC DSYMMatVec
#else
#define MYSYMMATVEC SSYMMatVec
#endif
// endif-else _DOUBLE_REAL_

#else
// !defined _DOUBLE_REAL_ and !defined _SINGLE_REAL_

#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)
#define TRANSA        "t"

#ifdef _SINGLE_COMPLEX_
#define MYSYMMATVEC CSYMMatVec
#else
#define MYSYMMATVEC ZSYMMatVec
#endif
// end if-else _SINGLE_COMPLEX_

#else
// !defined _COMPLEX_SYMMETRIC_

#define CONJG(A)     (-(A))
#define TRANSA        "c"

#ifdef _SINGLE_COMPLEX_
#define MYSYMMATVEC CHERMatVec
#else
#define MYSYMMATVEC ZHERMatVec
#endif
// end if-else _SINGLE_COMPLEX_


#endif
// end if-else !defined _COMPLEX_SYMMETRIC_

#endif
// end if-else defined _DOUBLE_REAL_ or defined _SINGLE_REAL_


void MYSYMMATVEC(SPARSEMATRIX *A, FLOAT *x, FLOAT *y) {

  integer i,j,k,n=A->nc;
  

  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
      y[i]=0.0;
#else
      y[i].r=y[i].i=0.0;
#endif
  }

  // multiplication with A
  for (i=0; i<n; i++) {
      for (j=0; j<A->ncol[i]; j++) {
	  k=A->rowind[i][j];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  y[k]+=A->val[i][j]*x[i];
#else
	  y[k].r+=A->val[i][j].r*x[i].r-A->val[i][j].i*x[i].i;
	  y[k].i+=A->val[i][j].r*x[i].i+A->val[i][j].i*x[i].r;
#endif
      }
  }
  // multiplication with CONJG((A-diag(A))^T)
  for (i=0; i<n; i++) {
      for (j=0; j<A->ncol[i]; j++) {
	  k=A->rowind[i][j];
	  if (k!=i) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     y[i]+=A->val[i][j]*x[k];
#else
	     y[i].r+=A->val[i][j].r*x[k].r-CONJG(A->val[i][j].i)*x[k].i;
	     y[i].i+=A->val[i][j].r*x[k].i+CONJG(A->val[i][j].i)*x[k].r;
#endif
	  }
      }
  }
  
} // end SYMMATVEC
