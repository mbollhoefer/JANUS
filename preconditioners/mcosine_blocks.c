/* $Id: mcosine_blocks.c 3620 2017-09-09 15:19:03Z bolle $ 

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


#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

void MCOSINE_BLOCKS(SPARSEMATRIX *A,
		    integer *p, integer *invq,
		    integer *blocksize, integer *nblocks,
		    REALS tau)
{  
  /*
    modified cosine-based compression
     
    Input
    -----
    A           sparse square matrix
    tau         tolerance

    Input/Output
    ------------
    p,invq      permutation vectors to group A(q,p) into blocks
                the given permutations are updated on exit
    Output
    ------
    blocksize   block size of the diagonal blocks in A(perm,perm)
    nblocks     number of blocks
  */


  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,tt,j_next,col, // counters
    
          *Ahead,      // head of the linked list for the rows of A
          *Alist,      // (incomplete) linked list for the rows of A
          *Afirst,     // physical position of the first nonzero entry of each
    
          *idxA,       // temporary index list of row indices
          cntA,        // length of idxA

          *group,      // check mark array for linking rows/column to a common block
          *count,      // counter for the formal scalar product
    
          *nnz_r,      // number of nonzeros for each row of A(q,p)
          *perm,       // permutation vector assoicated with blocks
          cnt;         // current length of perm
  REALS   expv,varr,varc,stddr,stddc;
  
  // index list 
  idxA=(integer *)malloc((size_t)n*sizeof(integer));

  // sort columns of A(q,p) such that for every column
  // of the permuted matrix the row indices and their
  // associated numerical values are taken in increasing
  // order
  // this is required to make the linked list arrays Ahead/Alist/Afirst
  // work properly when tracking the associated rows
  // Furthermore, this accelerates the cosine-based algorithm (and is required by the
  // implementation)
  // first switch to FORTRAN notation, because of the ILUPACK driver QQSORT
  for (j=0; j<n; j++)
      for (k=0;k<A->ncol[j];k++)
	  (A->rowind[j][k])++;
  for (j=0; j<n; j++)
      QQSORT(A->val[j],A->rowind[j],idxA,&(A->ncol[j]),invq);
  // switch back to C notation
  for (j=0; j<n; j++)
      for (k=0;k<A->ncol[j];k++)
	  (A->rowind[j][k])--;
#ifdef PRINT_INFO
  printf("matrix reordered\n");
  fflush(stdout);
#endif
  
  // set up arrays for A^T (A is stored by by columns)
  // implicitly these structures work with A(q,p)
  // head of the linked list for the rows of A
  Ahead =(integer *)malloc((size_t)n*sizeof(integer));
  // linked list for the leading rows of A up to step j
  Alist =(integer *)malloc((size_t)n*sizeof(integer));
  // position of the first entry at step j
  Afirst=(integer *)malloc((size_t)n*sizeof(integer));
  // clear head of the linked list (empty) 
  for (j=0; j<n; j++)
      Ahead[j]=-1;
  // init structures
  for (j=0; j<n; j++) {
      jj=p[j];
      // we have to make sure that we are still inside column p[j]
      if (0<A->ncol[jj]) {
	 // first row index k of A(q[k],p[j])
	 k=invq[A->rowind[jj][0]];
	 // position of the current first nonzero entry in column j
	 Afirst[j]=0;
	 // add new entry to the head of the list
         Alist[j]=Ahead[k];
	 Ahead[k]=j;
      } // end if
  } // end for j
  /*
  for (i=0; i<n; i++)
    printf("%8ld",Ahead[i]);
  printf("\n");
  for (i=0; i<n; i++)
    printf("%8ld",Alist[i]);
  printf("\n");
    for (i=0; i<n; i++)
    printf("%8ld",Afirst[i]);
  printf("\n");
  */
#ifdef PRINT_INFO
  printf("matrix linked list set up\n");
  fflush(stdout);
#endif

  // groups to be gathered into a single block
  group=(integer *)malloc((size_t)n*sizeof(integer));
  for (i=0; i<n; i++)
      group[i]=-1;
  // formal scalar product
  count=(integer *)calloc(n,sizeof(integer));
  // additional permutation to compress blocks
  perm=(integer *)malloc((size_t)n*sizeof(integer));
   
  // nonzeros per row in A(q,p)
  nnz_r=(integer *)calloc(n,sizeof(integer));
  for (j=0; j<n; j++) {
      jj=p[j]; // shortcut
      for (k=0; k<A->ncol[jj]; k++) {
	  // row count via accessing the columns quickly
	  ii=A->rowind[jj][k];
	  // q[i]=ii, i.e. we consider A(q[i],p[j])
	  i=invq[ii];
	  nnz_r[i]++;
      } // end for k	  
  } // end for j

  // compute arithmetic mean and standard deviations
  expv=varr=varc=0.0;
  for (j=0; j<n; j++) {
      expv+=nnz_r[j];
      varr+=nnz_r[j]*nnz_r[j];
      varc+=A->ncol[j]*A->ncol[j];
  } // end for j
  expv/=n;
  stddr=sqrt((varr-n*expv*expv)/n);
  stddc=sqrt((varc-n*expv*expv)/n);

  // exclude rows/columns with far too many nonzeros
  cnt=n-1;
  for (j=0; j<n; j++) {
      // check row/column j of A(q,p)
      jj=p[j];
      if (nnz_r[j]>=expv+2*stddr || A->ncol[jj]>=expv+2*stddc)
	 group[j]=++cnt;
  } // end for j

  
  
  
#ifdef PRINT_INFO
  printf("start main loop\n");
  fflush(stdout);
#endif
  
  // counter for the permutation vector
  cnt=0;
  //  counter for the current block
  *nblocks=-1;
  for (i=0; i<n; i++) {
      // if i-th row of A(q,p) is not yet assigned to any block
      if (group[i]<0) {            
	 // start a new block
	 blocksize[++(*nblocks)]=1;
	 // store permutation
	 perm[cnt++]=i;

	 // counter for the row index vector of relevant indices
	 cntA=0;
	 // find columns where the nonzeros A(q[i],p) are via linked list
	 ii=p[i]; // shortcut
	 // check ROW i of A(q,p)	    
	 // since A is stored by columns we access row i via the linked list
	 j=Ahead[i];
	 while (j>=0) {
	       // printf("scanning column %ld\n",j);
	       jj=p[j]; // shortcut

	       // ignore dense columns
	       if (group[j]<n) {
		  // scan sub-diagonal row indices of A(q,p[j]) starting from the back
		  for (k=A->ncol[jj]-1; k>=0; k--) {
		      // current index
		      l=A->rowind[jj][k];
		      // q[col]=l, i.e. we consider A(q[col],p[j])
		      col=invq[l];
		      // diagonal entry reached
		      if (col<=i)
			 break;
		      // otherwise, if col is not yet assigned
		      if (group[col]<0) {
			 // store new relevant index
			 if (count[col]==0) 
			    idxA[cntA++]=col;
			 // increase logical scalar product
			 count[col]++;
		      } // end if
		  } // end for k
	       } // end if group[j]<n
	       
	       // advance to the next column associated with row i
	       j=Alist[j];
	 } // end while j>=0
       
	 // scan again all relevant indices after their scalar product
	 // has been fully assembled
	 for (j=0; j<cntA; j++) {
	     col=idxA[j];
	     // consider block if the cosine of the angle is large enough
	     if (count[col]*count[col]>tau*tau*nnz_r[i]*nnz_r[col]) {
	        group[col]=i;
		// increase block size of the current block
		blocksize[*nblocks]++;
		// store permutation
		perm[cnt++]=col;
	     } // end if
	     // clear scalar product array
	     count[col]=0;
	 } // end for col
      } // end if group[i]<0
      
      // ---------------------------------------------------------------
      // ----- update linked list for the next row  greater than i -----
      // beginning of the linked list
      j=Ahead[i];
      // while list is non-empty
      while (j>=0) {
	    // position of the current leading entry in column j referring to A(q[i:end],p[j])
	    m=Afirst[j];
	    jj=p[j]; // shortcut
	    // row index tt of A(tt,p[j])
	    tt=A->rowind[jj][m];
	    // q[t]=tt,  we should have t=i
	    t=invq[tt];
	    // increment position of the next leading entry in column j
	    Afirst[j]=++m;
	    // buffer next column to be considered
	    j_next=Alist[j];
	    // we have to make sure that we are still inside column j
	    if (m<A->ncol[jj]) {
	       // row index tt of A(tt,p[j]) we should have t=i
	       tt=A->rowind[jj][m];
	       // q[t]=tt
	       t=invq[tt];
	       // add new entry to the head of the list
	       Alist[j]=Ahead[t];
	       Ahead[t]=j;
	    } // end if 
	    // restore next column
	    j=j_next;
      } // end while j>=0      
      // -----------------------------------------------
      // --- END update linked list for the next row ---
      /*
	for (m=0; m<n; m++)
	    printf("%8ld",Ahead[m]);
	printf("\n");
	for (m=0; m<n; m++)
	    printf("%8ld",Alist[m]);
	printf("\n");
	for (m=0; m<n; m++)
	    printf("%8ld",Afirst[m]);
	printf("\n");
      */

  } // end for i
#ifdef PRINT_INFO
  printf("main loop completed\n");
  fflush(stdout);
#endif

  // add rows/columns with far too many nonzeros at the end
  for (j=0; j<n; j++) {
      // check row/column j of A(q,p)
      if (group[j]>=n) {
	 ++(*nblocks);
	 // set block size of the current block to 1
	 blocksize[*nblocks]=1;
	 // store permutation
	 perm[cnt++]=j;
	 group[j]=-1;
      } // end if
  } // end for j
  // switch nblocks from index position to number of blocks
  ++(*nblocks);

  // re-arrange permutation p <- p[perm], q <- q[perm], invq <- invperm[invq]
  for (i=0; i<n; i++)
      count[i]=p[perm[i]]; // use count as buffer for p[perm]
  memcpy(p,count,(size_t)n*sizeof(integer));
  for (i=0; i<n; i++)
      group[perm[i]]=i; // use group as buffer for invperm
  for (i=0; i<n; i++)
      count[i]=group[invq[i]]; // use count as buffer for invperm[invq]
  memcpy(invq,count,(size_t)n*sizeof(integer));

  // sort original matrix back such that row indices and their
  // values are taken in increasing order
  for (j=0; j<n; j++)
      QSORT(A->val[j],A->rowind[j],idxA,&(A->ncol[j]));

#ifdef PRINT_INFO
  printf("matrix permuted back\n");
  fflush(stdout);
#endif

  free(Ahead);
  free(Alist);
  free(Afirst);
  free(group);
  free(count);
  free(perm);
  free(idxA);
  free(nnz_r);
} // cosine_blocks
