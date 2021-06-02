/* $Id: janussol.c 4161 2018-04-26 09:24:41Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <blas.h>

#include <janus.h>
#include <ilupackmacros.h>


void JanusSol(JanusPrec *PREC, void *x, void *y, void *buff, integer m) {

  integer i,j,k,n=PREC->n, *p=PREC->p, *invq=PREC->invq, dummy=-(PREC->n+1);
  double *SL=PREC->SL, *SR=PREC->SR, *dx, *dy, *dz;
  doublecomplex *zx, *zy, *zz;
  
  if (SR==NULL)
     SR=SL;


  if (PREC->isreal) {
     dx=(double *)x;
     dy=(double *)y;
     if (buff==NULL)
        dz=(double *)malloc((size_t)n*m*sizeof(double));
     else
        dz=buff;
     // left scaling + permutation
     if (SL!=NULL)
        for (j=0; j<m; j++,dx+=n,dz+=n)
	    for (i=0; i<n; i++) {
	        k=invq[i];  
		dz[k]=SL[i]*dx[i];
	    } // end for i
     else // SL==NULL
        for (j=0; j<m; j++,dx+=n,dz+=n)
	    for (i=0; i<n; i++) {
	        k=invq[i];  
		dz[k]=dx[i];
	    } // end for i
     dx-=m*n;
     dz-=m*n;
     
     if (PREC->issymmetric || PREC->ishermitian) {
        if (PREC->isdefinite)
	   if (PREC->invert_blocks) {
	      // printf("call symbilusol with inverse diagonal block\n");
	      DSYMbilusol((DSparseBlockMatrix *)(PREC->BL),
			  (DSparseBlockMatrix *)(PREC->BiD),
			  (DSparseBlockMatrix *)(PREC->BUT),
			  NULL, dz, dy, n, m);
	   }
	   else {
	      if (PREC->pivots==NULL) {
		 // printf("call symbilusol with factorized Cholesky diagonal block\n");
		 DSYMbilusol((DSparseBlockMatrix *)(PREC->BL),
			     (DSparseBlockMatrix *)(PREC->BiD),
			     (DSparseBlockMatrix *)(PREC->BUT),
			     &dummy, dz, dy, n, m);
	      }
	      else {
		 // printf("call symbilusol with factorized LDL^T diagonal block\n");
		 DSYMbilusol((DSparseBlockMatrix *)(PREC->BL),
			     (DSparseBlockMatrix *)(PREC->BiD),
			     (DSparseBlockMatrix *)(PREC->BUT),
			     PREC->pivots, dz, dy, n, m);
	      }
	   }
	else {
	   /*
	   if (PREC->pivots==NULL)
	      printf("call symbilusol with inverse diagonal block\n");
	   else
	      printf("call symbilusol with factorized LDL^T block\n");
	   */
	   DSYMbilusol((DSparseBlockMatrix *)(PREC->BL),
		       (DSparseBlockMatrix *)(PREC->BiD),
		       (DSparseBlockMatrix *)(PREC->BUT),
		       PREC->pivots, dz, dy, n, m);
	}
     }
     else {
        Dbilusol((DSparseBlockMatrix *)(PREC->BL),
		 (DSparseBlockMatrix *)(PREC->BiD),
		 (DSparseBlockMatrix *)(PREC->BUT),
		 PREC->pivots, dz, dy, n, m);
     }
     // inverse permutation + right scaling
     if (SR!=NULL)
        for (j=0; j<m; j++,dy+=n,dz+=n)
	    for (i=0; i<n; i++) {
	        k=p[i];  
		dy[k]=dz[i];
		dy[k]*=SR[k];
	    } // end for i
     else
        for (j=0; j<m; j++,dy+=n,dz+=n)
	    for (i=0; i<n; i++) {
	        k=p[i];  
		dy[k]=dz[i];
	    } // end for i
     dy-=m*n;
     dz-=m*n;
     
     if (buff==NULL)
        free(dz);
  }
  else { // complex-valued case
     zx=(doublecomplex *)x;
     zy=(doublecomplex *)y;
     if (buff==NULL)
        zz=(doublecomplex *)malloc(n*sizeof(doublecomplex));
     else
        zz=buff;
     // left scaling + permutation
     if (SL!=NULL)
        for (j=0; j<m; j++,zx+=n,zz+=n)
	    for (i=0; i<n; i++) {
	        k=invq[i];  
		zz[k].r=SL[i]*zx[i].r;
		zz[k].i=SL[i]*zx[i].i;
	    } // end for i
     else
        for (j=0; j<m; j++,zx+=n,zz+=n)
	    for (i=0; i<n; i++) {
	        k=invq[i];  
		zz[k].r=zx[i].r;
		zz[k].i=zx[i].i;
	    } // end for i
     zx-=m*n;
     zz-=m*n;
     if (PREC->issymmetric) {
        ZSYMbilusol((ZSparseBlockMatrix *)(PREC->BL),
		    (ZSparseBlockMatrix *)(PREC->BiD),
		    (ZSparseBlockMatrix *)(PREC->BUT),
		    PREC->pivots, zz, zy, n, m);
     }
     else if (PREC->ishermitian) {
        if (PREC->isdefinite)
	   if (PREC->invert_blocks)
	      ZHERbilusol((ZSparseBlockMatrix *)(PREC->BL),
			  (ZSparseBlockMatrix *)(PREC->BiD),
			  (ZSparseBlockMatrix *)(PREC->BUT),
			  NULL, zz, zy, n, m);
	   else
	      ZHERbilusol((ZSparseBlockMatrix *)(PREC->BL),
			  (ZSparseBlockMatrix *)(PREC->BiD),
			  (ZSparseBlockMatrix *)(PREC->BUT),
			  &dummy, zz, zy, n, m);
	else
	   ZHERbilusol((ZSparseBlockMatrix *)(PREC->BL),
		       (ZSparseBlockMatrix *)(PREC->BiD),
		       (ZSparseBlockMatrix *)(PREC->BUT),
		       PREC->pivots, zz, zy, n, m);
     }
     else {
        Zbilusol((ZSparseBlockMatrix *)(PREC->BL),
		 (ZSparseBlockMatrix *)(PREC->BiD),
		 (ZSparseBlockMatrix *)(PREC->BUT),
		 PREC->pivots, zz, zy, n, m);
     }
     // inverse permutation + right scaling
     if (SR!=NULL)
        for (j=0; j<m; j++,zy+=n,zz+=n)
	    for (i=0; i<n; i++) {
	        k=p[i];  
		zy[k]=zz[i];
		zy[k].r*=SR[k];
		zy[k].i*=SR[k];
	    } // end for i
     else // SR==NULL
        for (j=0; j<m; j++,zy+=n,zz+=n)
	    for (i=0; i<n; i++) {
	        k=p[i];  
		zy[k]=zz[i];
	    } // end for i
     zy-=m*n;
     zz-=m*n;
     if (buff==NULL)
        free(zz);
  } 
} // end JanusSol
