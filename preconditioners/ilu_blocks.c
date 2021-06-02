/* $Id: ilu_blocks.c 3664 2017-09-28 10:46:32Z bolle $ */
/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	September 28, 2017. JANUS Block ILU R1.0.  

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
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>




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

integer ILU_BLOCKS(SPARSEMATRIX *A,
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
          *Asubdiag,    // array of sub-diagonal entries
          *Asuperdiag;  // array of super-diagonal entries

  REALS   bnd,bndl,bndr, // bounds, parameters used for inverse of 2x2 blocks
          a11,a12,a21,a22;   // auxiliary variables for determinant
  

  

  // diagonal entries 
  Adiag     =(FLOAT *)malloc(n*sizeof(FLOAT));
  // sub-diagonal entries 
  Asubdiag  =(FLOAT *)malloc(n*sizeof(FLOAT));
  // super-diagonal entries 
  Asuperdiag=(FLOAT *)malloc(n*sizeof(FLOAT));
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
  for (i=0; i<n; i++)
      Adiag[i]=Asubdiag[i]=Asuperdiag[i]=0.0;
#else
  for (i=0; i<n; i++)
      Adiag[i].r=Asubdiag[i].r=Asuperdiag[i].r
	=Adiag[i].i=Asubdiag[i].i=Asuperdiag[i].i=0.0;
#endif



  // scan A(p,p) in order to find the diagonal entry and super/sub diagonal entry
  // though they exist
  if (SL==NULL && SR==NULL) {
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

	     // check for existing diagonal entry A(p[i],p[i])
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asubdiag[i]=val;
	     // check for the super-diagonal entry A(p[i-1],p[i])
	     else if (k==i-1)
	        Asuperdiag[i-1]=val;
	 } // end for j
     } // end for i
  }
  else if (SL==NULL) { // SR!=0
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
	     // SL(p[k],p[k]) A(p[k],p[i]) SR(p[i],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     val=pA[j]*SR[ii];
#else
	     val.r=pA[j].r*SR[ii];
	     val.i=pA[j].i*SR[ii];
#endif
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asubdiag[i]=val;
	     // check for the super-diagonal entry A(p[i-1],p[i])
	     else if (k==i-1)
	        Asuperdiag[i-1]=val;
	 } // end for j
     } // end for i
  } // end if-else-if SL=0
  else if (SR==NULL) { // SL!=0
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
	     // SL(p[k],p[k]) A(p[k],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     val=SL[kk]*pA[j];
#else
	     val.r=SL[kk]*pA[j].r;
	     val.i=SL[kk]*pA[j].i;
#endif
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asubdiag[i]=val;
	     // check for the super-diagonal entry A(p[i-1],p[i])
	     else if (k==i-1)
	        Asuperdiag[i-1]=val;
	 } // end for j
     } // end for i
  } // end if-else-if SR=0
  else { // SL!=0 and SR!=0
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
	     // SL(p[k],p[k]) A(p[k],p[i]) SR(p[i],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     val=SL[kk]*pA[j]*SR[ii];
#else
	     val.r=SL[kk]*pA[j].r*SR[ii];
	     val.i=SL[kk]*pA[j].i*SR[ii];
#endif
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asubdiag[i]=val;
	     // check for the super-diagonal entry A(p[i-1],p[i])
	     else if (k==i-1)
	        Asuperdiag[i-1]=val;
	 } // end for j
     } // end for i
  } // end if-else-if-else SL!=0 and SR!=0


 
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
	   // [a11 a13] = |A(p[j:j+1],p[j:j+1])|
	   // [a21 a22] 
	   // a11=|A(p[j],p[j])|
	   a11=bnd;
	   // a12=|A(p[j],p[j+1])|
	   a12=FABS(Asuperdiag[j]);
	   // a21=|A(p[j+1],p[j])|
	   a21=FABS(Asubdiag[j]);
	   // a22=|A(p[j+1],p[j+1])|
	   a22=FABS(Adiag[j+1]);
	   // check whether blocking columns j:j+1 is superior
	   // only take this into account, if the off-diagonal entries are
	   // sufficiently large
	   if (MAX(a12,a21)>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
	      // use 1/||A(p[j:j+1],p[j:j+1])^{-1}|| as measure
	      // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
	     val=Adiag[j]*Adiag[j+1]-Asuperdiag[j]*Asubdiag[j];
#else
	     val.r=(Adiag[j].r     *Adiag[j+1].r -Adiag[j].i     *Adiag[j+1].i)
	          -(Asuperdiag[j].r*Asubdiag[j].r-Asuperdiag[j].i*Asubdiag[j].i);
	     val.i=(Adiag[j].r     *Adiag[j+1].i +Adiag[j].i     *Adiag[j+1].r)
	          -(Asuperdiag[j].r*Asubdiag[j].i+Asuperdiag[j].i*Asubdiag[j].r);
#endif
	     // 1/||A(p[j:j+1],p[j:j+1])^{-1}||
	     bndr=FABS(val)/sqrt(MAX(a11+a21,a12+a22)*MAX(a11+a12,a21+a22));

	     // check if a 2x2 pivot is preferred, otherwise skip it
	     if (bndr>TWO_BY_TWO_BOUND*a11) {
	        // use column j+1 as additional column
	        colj=j+1;
	     } // end if
	   } // end if a12,a21 are large enough
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
  // printf("sub-diagonal part\n");
  // for (i=0; i<n-1; i++)
  //    printf("%12.4le",Asubdiag[i]);
  // printf("\n");
  // fflush(stdout);
  // printf("super-diagonal part\n");
  // for (i=0; i<n-1; i++)
  //    printf("%12.4le",Asuperdiag[i]);
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
  free(Asubdiag);
  free(Asuperdiag);

  
  return (n-j);
} // end ilu_blocks
