/* $Id: matvec.c 3608 2017-09-05 07:33:23Z bolle $ */
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

void MATVEC(SPARSEMATRIX *A, FLOAT *x, FLOAT *y) {

  integer i,j,k,n=A->nc,l,*idx;
  FLOAT *p, val;

  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
      y[i]=0.0;
#else
      y[i].r=y[i].i=0.0;
#endif
  }

  for (i=0; i<n; i++) {

      p=A->val[i];
      val=x[i];
      idx=A->rowind[i];
      l=A->ncol[i];
#ifdef _USE_MKL_
      // use sparse AXPYI y[A->rowind[i]]= x[i] * A->val[i] + y[A->rowind[i]]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      CBLAS_AXPYI(l, val,  p,idx, y); 
#else
      CBLAS_AXPYI(l, &val, p,idx, y); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
		    
#else // hand-coded loop
      for (j=0; j<l; j++) {
	  // position of k in y
	  k=idx[j];
	  // downdate buff - BL{i}*BUT{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	  y[k]+=(*p++)*val;
#else
	  y[k].r+=p->r*val.r-p->i*val.i;
	  y[k].i+=p->r*val.i+p->i*val.r;
	  p++;
#endif
      } // end for j
#endif //-else _USE_MKL_

      /*
      for (j=0; j<A->ncol[i]; j++) {
	  k=A->rowind[i][j];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  y[k]+=A->val[i][j]*x[i];
#else
	  y[k].r+=A->val[i][j].r*x[i].r-A->val[i][j].i*x[i].i;
	  y[k].i+=A->val[i][j].r*x[i].i+A->val[i][j].i*x[i].r;
#endif
      }
      */
  }
} // end MATVEC
