/* $Id: janusdelete.c 6264 2020-05-15 08:52:52Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

void JanusDelete(JanusPrec *PREC) {

   if (PREC==NULL)
      return;
   PREC->SL    =FREE(PREC->SL);
   PREC->SR    =FREE(PREC->SR);
   PREC->p     =FREE(PREC->p);
   PREC->invq  =FREE(PREC->invq);
   PREC->pivots=FREE(PREC->pivots);

   if (PREC->isreal) {
      if (PREC->BL!=NULL) {
	 DSparseBlockDelete((DSparseBlockMatrix *)PREC->BL);
	 PREC->BL =FREE(PREC->BL);
      }
      if (PREC->BiD!=NULL) {
	 DSparseBlockDelete((DSparseBlockMatrix *)PREC->BiD);
	 PREC->BiD=FREE(PREC->BiD);
      }
      if (PREC->BUT!=NULL) {
	 DSparseBlockDelete((DSparseBlockMatrix *)PREC->BUT);
	 PREC->BUT=FREE(PREC->BUT);
      }
   }
   else {
      if (PREC->BL!=NULL) {
	 ZSparseBlockDelete((ZSparseBlockMatrix *)PREC->BL);
	 PREC->BL =FREE(PREC->BL);
      }
      if (PREC->BiD!=NULL) {
	 ZSparseBlockDelete((ZSparseBlockMatrix *)PREC->BiD);
	 PREC->BiD=FREE(PREC->BiD);
      }
      if (PREC->BUT!=NULL) {
	 ZSparseBlockDelete((ZSparseBlockMatrix *)PREC->BUT);
	 PREC->BUT=FREE(PREC->BUT);
      }

   }
   PREC->BL=PREC->BiD=PREC->BUT=NULL;
     
   PREC->n=0;
   PREC->isdefinite=PREC->isreal=PREC->issymmetric=PREC->ishermitian=0;
}


void SparseDelete(SparseMatrix *A) {

   if (A==NULL)
      return;
   if (A->isreal) {
     /*
      if (A->issingle) 
         SSparseDelete((SSparseMatrix *)A);
      else
     */
         DSparseDelete((DSparseMatrix *)A);
   }
   else {
     /*
      if (A->issingle) 
	 CSparseDelete((CSparseMatrix *)A);
      else
     */
	 ZSparseDelete((ZSparseMatrix *)A);

   }
}

void SparseBlockDelete(SparseBlockMatrix *A) {

   if (A==NULL)
      return;
   if (A->isreal) {
     /*
      if (A->issingle) 
         SSparseBlockDelete((SSparseBlockMatrix *)A);
      else
     */
         DSparseBlockDelete((DSparseBlockMatrix *)A);
   }
   else {
     /*
      if (A->issingle) 
	 CSparseBlockDelete((CSparseBlockMatrix *)A);
      else
     */
	 ZSparseBlockDelete((ZSparseBlockMatrix *)A);

   }
}
