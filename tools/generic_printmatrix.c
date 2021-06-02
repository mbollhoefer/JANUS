/* $Id: generic_printmatrix.c 3608 2017-09-05 07:33:23Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>


void PrintMatrix(SparseMatrix *A) {
  if (A->isreal)
     DPrintMatrix((DSparseMatrix *)A);
  else
     ZPrintMatrix((ZSparseMatrix *)A);
} // end PrintMatrix
