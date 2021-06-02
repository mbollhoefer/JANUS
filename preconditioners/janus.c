/* $Id: janus.c 6318 2020-06-01 15:35:35Z bolle $ 

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	August 28, 2017. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <janus.h>
#include <ilupackmacros.h>




integer JanusFactor(SparseMatrix *A, JanusPrec *PREC, JanusOptions options)
{  
  integer ierr=0,n=A->nr;

  
  PREC->n=n;
  PREC->SL  =(void *)malloc(n*sizeof(double));
  
  PREC->p   =(integer *)malloc(n*sizeof(integer));
  PREC->invq=(integer *)malloc(n*sizeof(integer));
  PREC->invert_blocks=options.invert_blocks;
  if (options.invert_blocks==0 && !A->isdefinite)
     PREC->pivots=(integer *)malloc(n*sizeof(integer));
  else
     PREC->pivots=NULL;
  
  if (A->isreal) {
     
     PREC->BL =(SparseBlockMatrix *)malloc(1*sizeof(DSparseBlockMatrix));
     PREC->BiD=(SparseBlockMatrix *)malloc(1*sizeof(DSparseBlockMatrix));
     PREC->isreal=1;
     if (A->issymmetric || A->ishermitian) {
        PREC->BUT=NULL;
	PREC->SR=NULL;
	PREC->issymmetric=1;
	PREC->ishermitian=1;
	if (A->isdefinite) {
	   PREC->isdefinite=1;
	   ierr=DSPDbilu_driver((DSparseMatrix *)A,
				(DSparseBlockMatrix *) (PREC->BL),
				(DSparseBlockMatrix *) (PREC->BiD),
				(DSparseBlockMatrix *) (PREC->BUT),
				(double *) (PREC->SL), (double *) (PREC->SR),
				PREC->p, PREC->invq,
				&(PREC->logdet),
				options);
       }
       else {
	  PREC->isdefinite=0;
	  ierr=DSYMbilu_driver((DSparseMatrix *)A,
			       (DSparseBlockMatrix *) (PREC->BL),
			       (DSparseBlockMatrix *) (PREC->BiD),
			       (DSparseBlockMatrix *) (PREC->BUT),
			       PREC->pivots,
			       (double *) (PREC->SL), (double *) (PREC->SR),
			       PREC->p, PREC->invq,
				&(PREC->logdet),
				&(PREC->isdefinite),
			       options);
       }
    }
    else {
       PREC->BUT=(SparseBlockMatrix *)malloc(1*sizeof(DSparseBlockMatrix));
       PREC->SR =(void *)malloc(n*sizeof(double));
       PREC->isdefinite=0;
       PREC->issymmetric=0;
       PREC->ishermitian=0;
       ierr=DGNLbilu_driver((DSparseMatrix *)A,
			    (DSparseBlockMatrix *) (PREC->BL),
			    (DSparseBlockMatrix *) (PREC->BiD),
			    (DSparseBlockMatrix *) (PREC->BUT),
			    PREC->pivots,
			    (double *) (PREC->SL), (double *) (PREC->SR),
			    PREC->p, PREC->invq,
			    &(PREC->logdet),
			    options);
    }
  }
  else { // non-real case
     
     PREC->BL  =(SparseBlockMatrix *)malloc(1*sizeof(ZSparseBlockMatrix));
     PREC->BiD =(SparseBlockMatrix *)malloc(1*sizeof(ZSparseBlockMatrix));
     PREC->isreal=0;
     if (A->issymmetric) {
        PREC->BUT=NULL;
	PREC->SR =NULL;
	PREC->issymmetric=1;
	PREC->ishermitian=0;
	PREC->isdefinite=0;
	ierr=ZSYMbilu_driver((ZSparseMatrix *)A,
			     (ZSparseBlockMatrix *) (PREC->BL),
			     (ZSparseBlockMatrix *) (PREC->BiD),
			     (ZSparseBlockMatrix *) (PREC->BUT),
			     PREC->pivots,
			     (double *) (PREC->SL), (double *) (PREC->SR),
			     PREC->p, PREC->invq,
			     &(PREC->logdet),
			     &(PREC->isdefinite),
			     options);
    }
    else if (A->ishermitian) {
       PREC->BUT=NULL;
       PREC->SR =NULL;
       PREC->issymmetric=0;
       PREC->ishermitian=1;
       if (A->isdefinite) {
	  PREC->isdefinite=1;
	  ierr=ZHPDbilu_driver((ZSparseMatrix *)A,
			       (ZSparseBlockMatrix *) (PREC->BL),
			       (ZSparseBlockMatrix *) (PREC->BiD),
			       (ZSparseBlockMatrix *) (PREC->BUT),
			       (double *) (PREC->SL), (double *) (PREC->SR),
			       PREC->p, PREC->invq,
			       &(PREC->logdet),
			       options);
       }
       else {
	  PREC->isdefinite=0;
	  ierr=ZHERbilu_driver((ZSparseMatrix *)A,
			       (ZSparseBlockMatrix *) (PREC->BL),
			       (ZSparseBlockMatrix *) (PREC->BiD),
			       (ZSparseBlockMatrix *) (PREC->BUT),
			       PREC->pivots,
			       (double *) (PREC->SL), (double *) (PREC->SR),
			       PREC->p, PREC->invq,
			       &(PREC->logdet),
			       &(PREC->isdefinite),
			       options);
       }
    }
    else {
       PREC->BUT=(SparseBlockMatrix *)malloc(1*sizeof(ZSparseBlockMatrix));
       PREC->SR =(void *)malloc(n*sizeof(double));
       PREC->isdefinite=0;
       PREC->issymmetric=0;
       PREC->ishermitian=0;
       ierr=ZGNLbilu_driver((ZSparseMatrix *)A,
			    (ZSparseBlockMatrix *) (PREC->BL),
			    (ZSparseBlockMatrix *) (PREC->BiD),
			    (ZSparseBlockMatrix *) (PREC->BUT),
			    PREC->pivots,
			    (double *) (PREC->SL), (double *) (PREC->SR),
			    PREC->p, PREC->invq,
			    &(PREC->logdet),
			    options);
    }
  }


  // block ILU failed
  if (ierr) {
     PREC->SL=FREE(PREC->SL);
     PREC->SR=FREE(PREC->SR);
  
     PREC->p   =FREE(PREC->p);
     PREC->invq=FREE(PREC->invq);

     PREC->pivots=FREE(PREC->pivots);
  
     
     PREC->BL =FREE(PREC->BL);
     PREC->BiD=FREE(PREC->BiD);
     PREC->BUT=FREE(PREC->BUT);
  } // end if ierr
  
  return (ierr);       

}
