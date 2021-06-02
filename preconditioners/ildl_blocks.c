/* $Id: ildl_blocks.c 3620 2017-09-09 15:19:03Z bolle $ */
/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	December 25, 2016. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2016 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYILDL_BLOCKS     SSYMildl_blocks
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYILDL_BLOCKS     DSYMildl_blocks
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYILDL_BLOCKS     CSYMildl_blocks
#define CONJG(A)       (A)
#else // double complex
#define MYILDL_BLOCKS     ZSYMildl_blocks
#define CONJG(A)       (A)
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYILDL_BLOCKS     SSYMildl_blocks
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYILDL_BLOCKS     DSYMildl_blocks
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYILDL_BLOCKS     CHERildl_blocks
#define CONJG(A)       (-(A))
#else // double complex
#define MYILDL_BLOCKS     ZHERildl_blocks
#define CONJG(A)       (-(A))
#endif //-if-elif-else single-real

#endif //-else complex-symmetric



#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define ELBOW_BLOCK     1.3333
#define BLOCK_EXT       4
#define TWO_BY_TWO_THRESHOLD 0.1
#define TWO_BY_TWO_BOUND     1.5

// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

integer MYILDL_BLOCKS(SPARSEMATRIX *A,
		      REALS *SL, REALS *SR, integer *p, integer *invq,
		      integer *blocksize, integer *nblocks,
		      REALS droptol)
{  
  /*
    given a reordering and scaling determine approximately a block structure
    with 1x1 and 2x2 pivots based on the smaller norm of their inverses


    Input
    -----
    A           assumed to be a nonsingular sparse matrix 
    p,invq      permutation p and inverse permutation invq associated with the a priori
                permutation matrix P
		we further assume that the row indices of A(p,p) in each column k=p(i)
		are already sorted in increasing order. 
    SL          real diagonal scaling matrix
    SR          not referenced

    Output
    ------
    blocksize   array with the initial block sizes of each diagonal block
                of the scaled and reordered system.
    nblocks     number of blocks within the array "blocksize"
  */

  integer n=A->nc,      // total size
          i,j,k,l, ii,colj,kk, // counters
          *pi,*pi2,    // index array pointers
    
          startblock,  // beginning and end of the current diagonal
          endblock;    // block
    
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2,     // temporary numerical pointers
          *Adiag,       // array of diagonal entries 
          *Asdiag;      // array of sub-diagonal entries

  REALS   bnd,bndl,bndr, // bounds, parameters used for inverse of 2x2 blocks
          a11,a21,a22;   // auxiliary variables for determinant
  

  

  // diagonal entries 
  Adiag  =(FLOAT *)malloc(n*sizeof(FLOAT));
  // sub-diagonal entries 
  Asdiag =(FLOAT *)malloc(n*sizeof(FLOAT));
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
  for (i=0; i<n; i++)
      Adiag[i]=Asdiag[i]=0.0;
#else
  for (i=0; i<n; i++)
      Adiag[i].r=Asdiag[i].r=Adiag[i].i=Asdiag[i].i=0.0;
#endif



  // scan A(p,p) in order to find the diagonal entry and super/sub diagonal entry
  // though they exist
  if (SL==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // A(p[k],p[i])
	     val=pA[j];
	     
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asdiag[i]=val;

	     // val <- CONJ(val)
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	     val.i=CONJG(val.i);
#endif
	     
	     // check for the super-diagonal entry A(p[i-1],p[i])=CONJ(A(p[i],p[i-1])
	     if (k==i-1)
	        Asdiag[i-1]=val;
	     
	 } // end for j
     } // end for i
  }
  else { // SL!=0
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // SL(p[k],p[k]) A(p[k],p[i]) SL(p[i],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     val=SL[kk]*pA[j]*SL[ii];
#else
	     val.r=SL[kk]*pA[j].r*SL[ii];
	     val.i=SL[kk]*pA[j].i*SL[ii];
#endif
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asdiag[i]=val;
	     
	     // val <- CONJ(val)
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	     val.i=CONJG(val.i);
#endif
	     
	     // check for the super-diagonal entry A(p[i-1],p[i])=CONJ(A(p[i],p[i-1])
	     if (k==i-1)
	        Asdiag[i-1]=val;
	     
	 } // end for j
     } // end for i
  } // end if-else SL=0


 
  // **************************************************************************
  // *****                       main loop                                *****
  // **************************************************************************
  *nblocks=0; // number of diagonal blocks of A
  j=0;
  while (j<n) {

#ifdef PRINT_INFO
        printf("column j=%3ld\n",j);
	fflush(stdout);
#endif

	// scalar bound
	bnd=FABS(Adiag[j]);
	// so far no additional column (i.e. scalar approach)
	colj=j;
	// check whether a 1x1 pivot or a 2x2 pivot should be chosen
	// possibly second check the subsequent column j+1
	if (j<n-1) {
#ifdef PRINT_INFO
	   printf("also check column j+1=%3ld\n",j+1);
	   fflush(stdout);
#endif
	   // [a11 a21] = |A(p[j:j+1],p[j:j+1])|
	   // [a21 a22] 
	   a11=bnd;
	   a21=FABS(Asdiag[j]);
	   a22=FABS(Adiag[j+1]);
	   // check whether blocking columns j:j+1 is superior
	   // only take this into account, if the sub-diagonal entry is
	   // sufficiently large
	   if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
	      // use 1/||A(p[j:j+1],p[j:j+1])^{-1}|| as measure
	      // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
	     val=Adiag[j]*Adiag[j+1]-Asdiag[j]*Asdiag[j];
#else
	     val.r=(Adiag[j].r *Adiag[j+1].r      -Adiag[j].i *Adiag[j+1].i)
	          -(Asdiag[j].r*Asdiag[j].r       -Asdiag[j].i*CONJG(Asdiag[j].i));
	     val.i=(Adiag[j].r *Adiag[j+1].i      +Adiag[j].i *Adiag[j+1].r)
	          -(Asdiag[j].r*CONJG(Asdiag[j].i)+Asdiag[j].i*Asdiag[j].r);
#endif
	     // 1/||A(p[j:j+1],p[j:j+1])^{-1}||
	     bndr=FABS(val)/MAX(a11+a21,a21+a22);

	     // check if a 2x2 pivot is preferred, otherwise skip it
	     if (bndr>TWO_BY_TWO_BOUND*a11) {
	        // use column j+1 as additional column
	        colj=j+1;
	     } // end if
	   } // end if a21 large enough
	} // end if j<n-1

	blocksize[(*nblocks)++]=colj-j+1;
	j=colj+1;

  } // end while j
  // **************************************************************************
  // *****                     END main loop                              *****
  // **************************************************************************


  // post processing
  j=0;
  for (i=0; i<*nblocks; i++)
      j+=blocksize[i];
#ifdef PRINT_INFO
  // printf("subdiagonal part\n");
  // for (i=0; i<n-1; i++)
  //    printf("%12.4le",Asdiag[i]);
  // printf("\n");
  // fflush(stdout);
  if (SL!=NULL) {
     printf("scaling\n");
     for (i=0; i<n; i++)
         printf("%12.4le",SL[i]);
     printf("\n");
     fflush(stdout);
  }
  printf("permutation\n");
  for (i=0; i<n; i++)
      printf("%6ld",p[i]);
  printf("\n");
  fflush(stdout);
#endif
#ifdef PRINT_INFO
  printf("block partitioning, number of blocks: %ld, check sum %ld\n",*nblocks,j);
  for (i=0; i<*nblocks; i++)
      printf("%4ld",blocksize[i]);
  printf("\n");
  fflush(stdout);
#endif


  // release memory
  free(Adiag);
  free(Asdiag);

  
  return (n-j);
} // end ildl_blocks

 

 
