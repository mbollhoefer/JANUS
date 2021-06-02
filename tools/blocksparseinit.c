/* $Id: blocksparseinit.c 7315 2021-05-28 21:00:20Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

void SPARSEBLOCKINIT(SPARSEBLOCKMATRIX *B) {
    B->rowind=NULL;
    B->colind=NULL;
    B->valD  =NULL;
    B->valE  =NULL;
    B->rowind=NULL;
    B->colind=NULL;
    B->nblockcol=NULL;
    B->nblockrow=NULL;
    B->nr=B->nc=B->nnz=B->nblocks=0;
}
