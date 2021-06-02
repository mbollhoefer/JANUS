/* $Id: janusnnz.c 5181 2019-07-10 19:21:25Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

#define MAX(A,B)        (((A)>(B))?(A):(B))

void JanusNnz(JanusPrec *PREC, integer *nnz, integer *mxblock, doubleprecision *avgblock, doubleprecision *stddev)
{

  integer i,k,l,m, nblocks;
  doubleprecision mu2;
  
  SparseBlockMatrix *BL =PREC->BL;
  SparseBlockMatrix *BiD=PREC->BiD;
  SparseBlockMatrix *BUT=PREC->BUT;

  if (BUT==NULL)
     BUT=BL;

  nblocks=BiD->nblocks;
  *nnz=0;
  *mxblock=0;
  *avgblock=0.0;
  mu2=0.0;
  for (i=0; i<nblocks; i++) {
      k=BiD->nblockcol[i];
      l=BL->nblockrow[i];
      m=BUT->nblockrow[i];
      *nnz+=k*(k+l+m);

      *mxblock=MAX(*mxblock,k);
      *avgblock+=k;
      mu2+=k*k;
  } /* end for i */

  // average block size
  *avgblock=*avgblock/nblocks;
  // standard deviation
  *stddev=sqrt((mu2-nblocks**avgblock**avgblock)/MAX(1,nblocks-1));
}
