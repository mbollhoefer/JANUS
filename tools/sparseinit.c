/* $Id: sparseinit.c 6094 2020-03-28 16:17:50Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>


void SPARSEINIT(SPARSEMATRIX *B) {
    B->nr=B->nc=B->nnz=0;
    B->isdefinite=B->isreal=B->issymmetric=B->ishermitian=B->isskew=B->issingle=0;
    B->rowind=NULL;
    B->val=NULL;
    B->ncol=NULL;
}

// $Id: sparseinit.c 6094 2020-03-28 16:17:50Z bolle $

