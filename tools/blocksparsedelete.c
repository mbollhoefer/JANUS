/* $Id: bilu.c 2480 2016-09-10 00:28:25Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

void SPARSEBLOCKDELETE(SPARSEBLOCKMATRIX *B) {
    integer i, nblocks=B->nblocks;

    if (B==NULL)
       return;
    for (i=0; i<nblocks; i++) {
        if (B->rowind!=NULL)
	   if (B->rowind[i]!=NULL)
	      free(B->rowind[i]);
        if (B->colind!=NULL)
	   if (B->colind[i]!=NULL)
	      free(B->colind[i]);
        if (B->valD!=NULL)
	   if (B->valD[i]!=NULL)
	      free(B->valD[i]);
        if (B->valE!=NULL)
	   if (B->valE[i]!=NULL)
	      free(B->valE[i]);        
    }
    if (B->rowind!=NULL)
       free(B->rowind);
    if (B->colind!=NULL)
       free(B->colind);
    if (B->valD!=NULL)
       free(B->valD);
    if (B->valE!=NULL)
       free(B->valE);

    if (B->nblockcol!=NULL)
       free(B->nblockcol);
    if (B->nblockrow!=NULL)
       free(B->nblockrow);
    B->nr=B->nc=B->nnz=B->nblocks=0;
}
