/* $Id: mattvec.c 3608 2017-09-05 07:33:23Z bolle $ */
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

void MATTVEC(SPARSEMATRIX *A, FLOAT *x, FLOAT *y) {

  integer i,j,k,n=A->nc,l,*idx;
  FLOAT *p, val;

  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
      y[i]=0.0;
#else
      y[i].r=y[i].i=0.0;
#endif
  } // end for i

  for (i=0; i<n; i++) {

      p=A->val[i];
      idx=A->rowind[i];
      l=A->ncol[i];
      
#ifdef _USE_MKL_
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      val=CBLAS_DOTI(l, p, idx, x);
#else
      CBLAS_DOTUI(l, p, idx, x, &val);
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
	   
#else // hand-coded loop   
      // value to collect scalar product (0)
      val=y[i];
      for (j=0; j<l; j++) {
	  k=idx[j];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  val+=p[j]*x[k];
#else
	  val.r+=p[j].r*x[k].r-p[j].i*x[k].i;
	  val.i+=p[j].r*x[k].i+p[j].i*x[k].r;
#endif
      } // end for j
#endif
      
      y[i]=val;
  } // end for i
} // end MatTVec
