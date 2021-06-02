/* $Id: janusinit.c 6263 2020-05-15 07:04:30Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

void JanusInit(JanusPrec *PREC) {

   PREC->SL    =NULL;
   PREC->SR    =NULL;
   PREC->p     =NULL;
   PREC->invq  =NULL;
   PREC->pivots=NULL;

   PREC->BL=PREC->BiD=PREC->BUT=NULL;
     
   PREC->n=0;
   PREC->invert_blocks=0;
   PREC->logdet.r=PREC->logdet.i=0.0;
   PREC->isdefinite=PREC->isreal=PREC->issymmetric=PREC->ishermitian=0;
}


void SparseInit(SparseMatrix *B) {
    B->nr=B->nc=B->nnz=0;
    B->isdefinite=B->isreal=B->issymmetric=B->ishermitian=B->isskew=B->issingle=0;
    B->rowind=NULL;
    B->val=NULL;
    B->ncol=NULL;
}

void SparseBlockInit(SparseBlockMatrix *B) {
    B->rowind=NULL;
    B->colind=NULL;
    B->valD  =NULL;
    B->valE  =NULL;
    B->rowind=NULL;
    B->colind=NULL;
    B->nblockcol=NULL;
    B->nblockrow=NULL;
    B->nr=B->nc=B->nblocks=0;
}
