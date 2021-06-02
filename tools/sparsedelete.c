/* $Id: sparsedelete.c 6264 2020-05-15 08:52:52Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>


void SPARSEDELETE(SPARSEMATRIX *B) {
    integer i, n;

    if (B==NULL)
       return;
    n=B->nc;
    for (i=0; i<n; i++) {
        if (B->rowind!=NULL)
	   if (B->rowind[i]!=NULL)
	      free(B->rowind[i]);
        if (B->val!=NULL)
	   if (B->val[i]!=NULL)
	      free(B->val[i]);
    }
    if (B->rowind!=NULL) {
       free(B->rowind);
       B->rowind=NULL;
    }
    if (B->val!=NULL) {
       free(B->val);
       B->val=NULL;
    }

    if (B->ncol!=NULL) {
       free(B->ncol);
       B->ncol=NULL;
    }
       
    B->nr=B->nc=B->nnz=0;
    B->isdefinite=B->isreal=B->issymmetric=B->ishermitian=B->isskew=B->issingle=0;
}

// $Id: sparsedelete.c 6264 2020-05-15 08:52:52Z bolle $

