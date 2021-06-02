/* $Id: generic_matvec.c 3608 2017-09-05 07:33:23Z bolle $ */
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

void MatVec(SparseMatrix *A, void *x, void *y) {

  if (A->isreal) {
     if (A->issymmetric || A->ishermitian) {
        DSYMMatVec((DSparseMatrix *)A, (double *)x, (double *)y);
     }
     else 
        DMatVec((DSparseMatrix *)A, (double *)x, (double *)y);
  }
  else {
     if (A->issymmetric) {
        ZSYMMatVec((ZSparseMatrix *)A, (doublecomplex *)x, (doublecomplex *)y);
     }
     else if (A->ishermitian) {
        ZHERMatVec((ZSparseMatrix *)A, (doublecomplex *)x, (doublecomplex *)y);
     }
     else 
        ZMatVec((ZSparseMatrix *)A, (doublecomplex *)x, (doublecomplex *)y);
  }
} // end MatVec
