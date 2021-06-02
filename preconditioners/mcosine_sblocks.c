/* $Id: mcosine_sblocks.c 3608 2017-09-05 07:33:23Z bolle $ */
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

void MCOSINE_SBLOCKS(SPARSEMATRIX *A,
		     integer *p, integer *invq,
		     integer *blocksize, integer *nblocks,
		     REALS tau)
{  
  /*
    modified cosine-based compression for symmetrically structured matrices
     
    Input
    -----
    A           sparse square symmetric matrix, only half of it being stored
    tau         tolerance

    Input/Output
    ------------
    p,invq      permutation vectors to group A(p,p) into blocks
                the given permutations are updated on exit
    Output
    ------
    blocksize   block size of the diagonal blocks in A(p,p)
    nblocks     number of blocks
  */


  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,tt,j_next,col, // counters
    
          *idxA,       // temporary index list of row indices
          cntA,        // length of idxA

          *group,      // check mark array for linking rows/column to a common block
          *count,      // counter for the formal scalar product

          *rowind, *rowindt,
    
          *nnz,*nnzt,        // number of nonzeros for each row/column of A(p,p)
          *perm,       // permutation vector assoicated with blocks
          cnt;         // current length of perm
        
  REALS   expv,var,stdd;
  SPARSEMATRIX AT;
  
  // index list 
  idxA=(integer *)malloc((size_t)n*sizeof(integer));

  // sort columns of A(p,p) such that for every column
  // of the permuted matrix the row indices and their
  // associated numerical values are taken in increasing
  // order
  // this is required to make the linked list arrays Ahead/Alist/Afirst
  // work properly when tracking the associated rows
  // Furthermore, this accelerates the cosine-based algorithm (and is required by the
  // implementation)
  // first switch to FORTRAN notation, because of the ILUPACK driver QQSORT
  nnz=A->ncol; // shortcut
  for (j=0; j<n; j++) {
    rowind=A->rowind[j]; // shortcut
      for (k=0;k<nnz[j];k++)
	  (rowind[k])++;
  } // end for j
  for (j=0; j<n; j++)
      QQSORT(A->val[j],A->rowind[j],idxA,&(nnz[j]),invq);
  // switch back to C notation
  for (j=0; j<n; j++) {
      rowind=A->rowind[j]; // shortcut
      for (k=0;k<nnz[j];k++)
	  (rowind[k])--;
  } // end for j
#ifdef PRINT_INFO
  printf("matrix reordered\n");
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
   
  // Compute pattern of (A(p,p)-diag(A(p,p)))^T, recall that only half of the 
  // symmetrically structured matrix is stored!
  nnzt=(integer *)calloc(n,sizeof(integer));
  AT.ncol=nnzt;
  AT.rowind=(integer **)malloc(n*sizeof(integer *));
  // first pass, count nonzeros for each row A(p(i),p) indirectly via columns
  for (j=0; j<n; j++) {
      // scan column p[j] of A
      jj=p[j]; // shortcut
      rowind=A->rowind[jj]; // shortcut
      
      for (k=0; k<nnz[jj]; k++) {
	  // row count via accessing the columns quickly
	  ii=rowind[k];
	  // p[i]=ii, i.e. we consider A(p[i],p[j])
	  i=invq[ii];
	  // exclude the diagonal part if it exists
	  if (i!=j)
	     nnzt[i]++;
      } // end for k	  
  } // end for j
  // allocate row index structure for the column of AT=A(p,p)^T
  for (j=0; j<n; j++) {
      AT.rowind[j]=(integer *)malloc(nnzt[j]*sizeof(integer));
      // reset nonzeros to re-use thema as counters
      nnzt[j]=0;
  } // end for j
  // second pass, insert indices
  for (j=0; j<n; j++) {
      // scan column p[j] of A
      jj=p[j]; // shortcut
      rowind=A->rowind[jj]; // shortcut
      for (k=0; k<nnz[jj]; k++) {
	  // row count via accessing the columns quickly
	  ii=rowind[k];
	  // p[i]=ii, i.e. we consider A(p[i],p[j])
	  i=invq[ii];
	  // exclude the diagonal part if it exists
	  if (i!=j) {
	     // counter for column i of AT
	     m=nnzt[i];
	     // note that this way, the indices in AT(:,i) are sorted in 
	     // increasing order automatically
	     AT.rowind[i][m]=j;
	     // increase counter, at the end it will denote the number of indices
	     nnzt[i]=m+1;
	  } // end if
      } // end for k	  
  } // end for j

  
  // compute arithmetic mean and standard deviation
  expv=var=0.0;
  for (j=0; j<n; j++) {
      jj=p[j]; // shortcut
      expv+=nnzt[j]+nnz[jj];
      var+=(nnzt[j]+nnz[jj])*(nnzt[j]+nnz[jj]);
  } // end for j
  expv/=n;
  stdd=sqrt((var-n*expv*expv)/n);

  // exclude rows/columns with far too many nonzeros
  cnt=n-1;
  for (j=0; j<n; j++) {
      // check row/column j of A(p,p)
      jj=p[j]; // shortcut
      if (nnzt[j]+nnz[jj]>=expv+2*stdd)
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
      // if i-th row of A(p,p) is not yet assigned to any block
      if (group[i]<0) {            
	 // start a new block
	 blocksize[++(*nblocks)]=1;
	 // store permutation
	 perm[cnt++]=i;

	 // counter for the row index vector of relevant indices
	 cntA=0;
	 // find columns where the nonzeros A(p[i],p) are
	 ii=p[i]; // shortcut
	 // check ROW i of A(p,p), i.e., check column AT(:,i)
	 for (m=0; m<nnzt[i]; m++) {
	     j=AT.rowind[i][m];
	     // printf("scanning column %ld\n",j);
	     jj=p[j]; // shortcut
	     rowind =A->rowind[jj]; // shortcut
	     rowindt=AT.rowind[j]; // shortcut

	     // ignore dense columns
	     if (group[j]<n) {
	        // scan sub-diagonal row indices of A(p,p[j]) starting from the back
	        for (k=nnz[jj]-1; k>=0; k--) {
		    // current index
		    l=rowind[k];
		    // p[col]=l, i.e. we consider A(p[col],p[j])
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
	        // scan sub-diagonal row indices of AT(:,j) starting from the back
	        for (k=nnzt[j]-1; k>=0; k--) {
		    // current index col of AT(col,j)
		    col=rowindt[k];
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
	 } // end for m

	 // check column i of A(p,p)
	 for (m=0; m<nnz[ii]; m++) {
	     jj=A->rowind[ii][m];
	     // p[j]=jj
	     j=invq[jj]; // shortcut
	     // printf("scanning column %ld\n",j);
	     rowind =A->rowind[jj]; // shortcut
	     rowindt=AT.rowind[j]; // shortcut

	     // ignore dense columns
	     if (group[j]<n) {
	        // scan sub-diagonal row indices of A(p,p[j]) starting from the back
	        for (k=nnz[jj]-1; k>=0; k--) {
		    // current index
		    l=rowind[k];
		    // p[col]=l, i.e. we consider A(p[col],p[j])
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
	        // scan sub-diagonal row indices of AT(:,j) starting from the back
	        for (k=nnzt[j]-1; k>=0; k--) {
		    // current index col of AT(col,j)
		    col=rowindt[k];
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
	 } // end for m
       
	 // scan again all relevant indices after their scalar product
	 // has been fully assembled
	 for (j=0; j<cntA; j++) {
	     col=idxA[j];
	     // consider block if the cosine of the angle is large enough
	     if (count[col]*count[col]>tau*tau*(nnz[i]+nnzt[i])*(nnzt[col]+nnz[col])) {
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

  } // end for i
  for (j=0; j<n; j++)
      free(AT.rowind[j]);
  free(AT.rowind);
  free(nnzt);
#ifdef PRINT_INFO
  printf("main loop completed\n");
  fflush(stdout);
#endif

  // add rows/columns with far too many nonzeros at the end
  for (j=0; j<n; j++) {
      // check row/column j of A(p,p)
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

  // re-arrange permutation p <- p[perm], invp <- invperm[invp]
  for (i=0; i<n; i++)
      count[i]=p[perm[i]]; // use count as buffer for p[perm]
  memcpy(p,count,(size_t)n*sizeof(integer));
  for (i=0; i<n; i++)
      group[perm[i]]=i; // use group as buffer for invperm
  for (i=0; i<n; i++)
      count[i]=group[invq[i]]; // use count as buffer for invperm[invp]
  memcpy(invq,count,(size_t)n*sizeof(integer));

  // sort original matrix back such that row indices and their
  // values are taken in increasing order
  for (j=0; j<n; j++)
      QSORT(A->val[j],A->rowind[j],idxA,&(nnz[j]));

#ifdef PRINT_INFO
  printf("matrix permuted back\n");
  fflush(stdout);
#endif

  free(group);
  free(count);
  free(perm);
  free(idxA);
} // mcosine_sblocks
