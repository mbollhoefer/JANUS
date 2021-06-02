/* $Id: bilu.c 7315 2021-05-28 21:00:20Z bolle $

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	September 29, 2017. JANUS Block ILU R1.0.  

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

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

// relaxation parameters when re-allocating buffers
#define ELBOW_BUFF      1.2
#define BUFF_EXT        100

// relaxation parameter when deciding whether to merge blocks or not
#define ELBOW_BLOCK     1.2
#define BLOCK_EXT       2

// upper bound for small size blocks s.t. level-3-BLAS do not pay off
#define SMALL_BLOCK     40.0

#define ABS_THRESHOLD   1e-2
#define REL_THRESHOLD   1e-1

#define PROGRESSIVE_AGGREGATION 1

// #define PRINT_INFO
// #define PRINT_CHECK
// #define printf mexPrintf 
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

integer BILU(SPARSEMATRIX *A,
	     SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT,
	     REALS *SL, REALS *SR, integer *p, integer *invq,
	     doublecomplex *determinant,
	     integer *blocksize, integer nblocks,
	     REALS droptol, integer ilu1t, integer pa, integer perturbation, 
	     integer symmetric_structure, integer *pivots, integer level_of_fill)
{  
  /*
    compute block incomplete LDU factorization with relative drop tolerance

    This routine approximately computes Q^T S_L A S_R P ~  L D U,
    where P, Q are given permutation matrices, S_L and S_R are real 
    diagnal scaling matrices, 
    L is lower unit block triangular matrix, 
    D is a block diagonal matrix and U is an upper unit block triangular
    matrix.

    Input
    -----
    A           assumed to be a nonsingular sparse matrix 
    p,invq      permutation p and inverse permutation invq associated with the a priori
                permutation matrices P and Q
    SL,SR       real diagonal scaling matrices
    blocksize   array with the initial block sizes of each diagonal block
                of the scaled and reordered system.
    nblocks     length of array blocksize		

    Output
    ------
    BL          lower block triangular matrix L with unit diagonal
    BiD         block diagonal matrix D^{-1} (LAPACK-factorized and inverted)
    BUT         transposed upper block triangular matrix U^T with unit diagonal
  */

  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,mm,tt,it,it_next,jt,jjt,jt_next,i2, // counters
          *invblock=NULL,// mapping scalar index j -> block i of A
          nrAL,ncAL,nrncAL,nrAU,ncAU,nrncAU,ncAD, // estimates for the blocks of A
          *pi,*pi2,    // index array pointers
          nbuff,nbuffC,// logical size of buff, buffC, buffDGLBLCK
          nbuffDGLBLCK,// 
          lda,ldb,ldc, // leading dimension used in combination with level-3-BLAS
          flag,        // flag indicating when a block column is exceeded
          bi,          // current block size
          nrowind,     // length of array Arowind
          ncolind,     // length of array Acolind
          i_bldu,      // index counter for the block ILU
          cnt_fill,    // counter for additional fill when merging blocks
          *rowind,ncol,*colind,nrow,// short cuts
    
          *idxL=NULL,   // temporary index list of row indices
          *idxposL=NULL,// associated inverse mapping storing the position
          cntL,        // length of idxL
          cntL_old,    // old length of idxL when needed

          *idxU=NULL,  // temporary index list of column indices
          *idxposU=NULL,// associated inverse mapping storing the position
          cntU,        // length of idxU
          cntU_old,    // old length of idxU when needed

          *idxbuff=NULL,// auxiliary index buffer
          skip_gthr_sctr, // flag when gather/scatter operation is not needed

          *Ahead=NULL, // head of the linked list for the rows of A
          *Alist=NULL, // (incomplete) linked list for the rows of A
          *Afirst=NULL,// physical position of the first nonzero entry of each 
                       // column as far as we still need to access this column
                       // to extract the associated rows
    
          *BLhead=NULL,// head of the linked list of block columns of BL
                       // required to extract the scalar rows of BL
          *BLlist=NULL,// (incomplete) linked list of block columns of BL
                       // required to extract the scalar rows of BL
          *BLfirst=NULL,// physical position of the first nonzero entry in each
                       // block column of BL as far as it is required to access
                       // a scalar row of BL
    
          *BUThead=NULL,// head of the linked list of block columns of BUT
                       // required to extract the scalar rows of BUT
          *BUTlist=NULL,// (incomplete) linked list of block columns of BUT
                       // required to extract the scalar rows of BUT
          *BUTfirst=NULL,// physical position of the first nonzero entry in each
                       // block column of BUT as far as it is required to access
                       // a scalar row of BUT
    
          *Arowind=NULL,// array with the nonzero row indices of the current
                       // sub-diagonal block of A
          *Acolind=NULL,// array with the nonzero column indices of the current
                       // transposed super-diagonal block of A

          flag_ilu1t=0,// flag for using ILU(1,droptol)
          mynblocks,   // copies of nblocks and blocksize, when ILU(1,droptol)
          *myblocksize,// is called, copied back on return
    
          ierr=0;      // return value bilu

  
  FLOAT   *AvalD=NULL, // buffer for the diagonal blocks of A
          *AvalL=NULL, // buffer for the current sub-diagonal block of A
          *AvalU=NULL,// buffer for the current transposed super-diagonal block of A
          *pL,*pD,*pU,*pb,*pC, // pointers to lower/diagonal/upper/buff/buffC
          val,v,ajj,    // temporary scalar numerical value
          *buff=NULL,  // buffer for the current block column/row of BL/BUT
          *buffC=NULL, // buffer for level-3-BLAS update
          *buffDGLBLCK=NULL,// buffer to hold a copy of the diagonal block in case
                       // the block LU decomposition needs to be re-computed, i.e.,
                       // when perturbation is turned on
          alpha, beta; // parameters for level-3-BLAS
  doublecomplex locdet;

  REALS   rval,rv,     // temporary scalar real numerical value
          as,rs,       // absolute and relative shift
          macheps=sqrt(GETEPS()),//square root of the machine precision 
          Anrm,        // infinity norm of the diagonal block
          Ajmax,       // maximum entry of |A(j+1:n,j)|
          Amax=0.0;    // maximum entry of |A|
  
  char    *transa, *transb; // strings for level-3-BLAS
  integer invert_blocks=1;

#ifdef _PROFILING_
  double timeBegin, timeBegin2,
         time_bilu=0.0,
         time_init_matrix=0.0,
         time_ilu1t=0.0,
         time_extract_block=0.0,
         time_LU_update_pass_1=0.0,
         time_LU_update_pass_2=0.0,
         time_scalar_update=0.0,
         time_rank_1_update=0.0,
         time_small_block_update=0.0,
         time_level_3_blas_update=0.0,
         time_level_3_blas_gather=0.0,
         time_level_3_blas_scatter=0.0,
         time_list_update=0.0,
         time_diagonal_block=0.0,
         time_off_diagonal_block=0.0,
         time_progressive_aggregation=0.0,
         time_post_processing=0.0;

  time_bilu=omp_get_wtime();
#endif

  determinant->r=0.0;
  determinant->i=0.0;

  // direct solver case
  if (droptol<=0.0)
     droptol=0.0;
  if (pivots!=NULL)
     invert_blocks=0;

  
  // compute max|A|
    /*
     for (j=0; j<n; j++) {
         jj=p[j];
         for (i=0; i<A->ncol[jj]; i++) 
	   printf("%8ld",invq[A->rowind[jj][i]]);
	 printf("\n");
         for (i=0; i<A->ncol[jj]; i++)
	   if (SL!=NULL && SR!=NULL)
	     printf("%8.1e",SL[A->rowind[jj][i]]*A->val[jj][i]*SR[jj]);
	   else
	     printf("%8.1e",A->val[jj][i]);
	 printf("\n");
     }
    */
     k=1;
     for (j=0; j<n; j++) {
         i=I_AMAX(&(A->ncol[j]),A->val[j],&k)-1;
	 l=A->rowind[j][i];
	 rval=FABS(A->val[j][i]);
	 // left scaling available
	 if (SL!=NULL)
	    rval*=SL[l];
	 // right scaling available
	 if (SR!=NULL)
	    rval*=SR[j];
	 Amax=MAX(Amax,rval);
     } // end for j
     // printf("Amax=%8.1le\n",Amax);
  
  // create auxiliary index buffer
  idxbuff=(integer *)malloc((size_t)n*sizeof(integer));
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc((size_t)(n+1)*sizeof(integer));
  idxL[0]=0;idxL++;
  // init with zeros
  idxposL=(integer *)calloc(n,sizeof(integer));
  idxU   =(integer *)malloc((size_t)(n+1)*sizeof(integer));
  idxU[0]=0;idxU++;
  // init with zeros
  idxposU=(integer *)calloc(n,sizeof(integer));

  // buffer for the current block column of BL/BUT
  nbuff =n;
  buff  =(FLOAT *)malloc((size_t)nbuff *sizeof(FLOAT));
  // buffer for level-3-BLAS
  nbuffC=n;
  buffC =(FLOAT *)malloc((size_t)nbuffC*sizeof(FLOAT));
  // buffer to hold a copy of the diagonal block if needed
  if (perturbation) {
     nbuffDGLBLCK=n;
     buffDGLBLCK=(FLOAT *)malloc((size_t)nbuffDGLBLCK*sizeof(FLOAT));
  }
  else {
     nbuffDGLBLCK=0;
     buffDGLBLCK=NULL;
  }

  

  // sort columns of A(q,p) such that for every column
  // of the permuted matrix the row indices and their
  // associated numerical values are taken in increasing
  // order
  // this is required to make the linked list arrays Ahead/Alist/Afirst
  // work properly when tracking the associated rows
  // switch to FORTRAN notation, because of the ILUPACK driver QQSORT
#ifdef _PROFILING_
  timeBegin = omp_get_wtime();
#endif
  for (j=0; j<n; j++) {
      rowind=A->rowind[j];
      ncol=A->ncol[j];
      for (k=0;k<ncol;k++)
	  (rowind[k])++;
  } // end for j
  for (j=0; j<n; j++)
      QQSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]),invq);
  // switch back to C notation
  for (j=0; j<n; j++) {
      rowind=A->rowind[j];
      ncol=A->ncol[j];
      for (k=0;k<ncol;k++)
	  (rowind[k])--;
  } // end for j
  
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
#ifdef _PROFILING_
  time_init_matrix=omp_get_wtime()-timeBegin;
#endif
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

  // no a priori block partitioning given, use ILU(1,droptol)
#ifdef _PROFILING_
  timeBegin = omp_get_wtime();
#endif
  if (nblocks==0 || blocksize==NULL) {
     flag_ilu1t=-1;
     myblocksize=blocksize;
     mynblocks=nblocks;
     blocksize=(integer *)malloc((size_t)n*sizeof(integer));
     if (ilu1t) {
        if (ilu1t==BLOCK_ILU1T) {
#ifdef PRINT_INFO
	   printf("call ILU1T\n");fflush(stdout);
#endif
	   if (symmetric_structure)
	      ierr=BILU1T_BLOCKS(A, SL,SR, p,invq, blocksize,&nblocks, droptol);
	   else
	      ierr=ILU1T_BLOCKS(A, SL,SR, p,invq, blocksize,&nblocks, droptol);
	}
	else if (ilu1t==BLOCK_ILUPT) {
#ifdef PRINT_INFO
	   printf("call ILUPT\n");fflush(stdout);
#endif
	   ierr=ILUPT_BLOCKS(A, SL,SR, p,invq, blocksize,&nblocks, droptol,level_of_fill);
	}
	else if (ilu1t==BLOCK_SUPERNODES) {
#ifdef PRINT_INFO
	   printf("call Supernode partitioning\n");fflush(stdout);
#endif
	   ierr=SUPERNODES(A, SL,SL, p,invq, blocksize,&nblocks, droptol);
	}
     }
     else {
#ifdef PRINT_INFO
        printf("ILU1T not required\n");fflush(stdout);
#endif
	if (symmetric_structure)
	   ierr=ILU_BLOCKS(A, SL,SR, p,invq, blocksize, &nblocks, droptol);
	else {
	   nblocks=n; 
	   for (i=0;i<n;i++)
	       blocksize[i]=1;
	   ierr=0;
	} // end if-else symmetric_structure
     }
     if (ierr) {
        // sort original matrix back such that row indices and their
        // values are taken in increasing order
        for (j=0; j<n; j++)
	    QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));

	idxL--;idxU=FREE(idxL);
	idxU--;idxU=FREE(idxU);
	idxbuff  =FREE(idxbuff  );
        idxposL  =FREE(idxposL  );
	idxposU  =FREE(idxposU  );
	buff     =FREE(buff     );
	buffC    =FREE(buffC    );
        Ahead    =FREE(Ahead    );
        Alist    =FREE(Alist    );
        Afirst   =FREE(Afirst   );
        blocksize=FREE(blocksize);
	
	if (perturbation)
	   buffDGLBLCK=FREE(buffDGLBLCK);

        return (ierr);
     }
  } // end if nblocks=0 or blocksize=0
#ifdef PRINT_INFO
  else {
     printf("ILU1T not required\n");fflush(stdout);
  }
#endif
#ifdef _PROFILING_
  time_ilu1t=omp_get_wtime()-timeBegin;
#endif

  
#ifdef _PROFILING_
  timeBegin=omp_get_wtime();
#endif
  // inverse mapping w.r.t. blocksize
  invblock=(integer *)malloc((size_t)n*sizeof(integer));
  j=0;
  for (i=0; i<nblocks; i++) {
      for (k=0; k<blocksize[i]; k++) {
  	  // scalar column j belongs to block column i
	  invblock[j++]=i;
      } // end for k
  } // end for i
  
  // set up sparse block data structures for BL,BUT,BiD
  BL->nc=BL->nr=n;
  BL->nblocks =nblocks;
  BL->nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BL->nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BL->colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BL->rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BL->valD     =NULL;
  BL->valE     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));

  BUT->nc=BUT->nr=n;
  BUT->nblocks =nblocks;
  BUT->nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BUT->nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BUT->colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BUT->rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BUT->valD     =NULL;
  BUT->valE     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));

  BiD->nc=BiD->nr=n;
  BiD->nblocks =nblocks;
  BiD->nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BiD->nblockrow=NULL;
  BiD->colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BiD->rowind   =NULL;
  BiD->valD     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));
  BiD->valE     =NULL;
  for (i=0; i<nblocks; i++) {
      BL->nblockrow[i]=BUT->nblockrow[i]=0;
      BL->nblockcol[i]=BUT->nblockcol[i]=0;
      BL->colind[i]=BL->rowind[i]=BUT->colind[i]=BUT->rowind[i]=BiD->colind[i]=NULL;
      BL->valE[i]=BUT->valE[i]=BiD->valD[i]=NULL;
  }
  // set up arrays for the rows BL^T (BL is stored by block columns)
  // these are required to access BL by rows during the computation
  // of the current block column of BUT
  // head of the linked list for the rows of BL
  BLhead =(integer *)malloc((size_t)nblocks*sizeof(integer));
  // linked list for the leading rows of BL up to step j
  BLlist =(integer *)calloc(nblocks,sizeof(integer));
  // position of the first entry at step j
  BLfirst=(integer *)calloc(nblocks,sizeof(integer));
  // clear head of the linked list (empty) 
  for (j=0; j<nblocks; j++)
      BLhead[j]=-1;

  // set up arrays for the rows BUT^T (BUT is stored by block columns)
  // these are required to access BUT by rows during the computation
  // of the current block column of BL
  // head of the linked list for the rows of BUT
  BUThead =(integer *)malloc((size_t)nblocks*sizeof(integer));
  // linked list for the leading rows of A up to step j
  BUTlist =(integer *)calloc(nblocks,sizeof(integer));
  // position of the first entry at step j
  BUTfirst=(integer *)calloc(n,sizeof(integer));
  // clear head of the linked list (empty) 
  for (j=0; j<nblocks; j++)
      BUThead[j]=-1;
  

  
  // buffer for the current diagonal block, block column/block row of A(q,p)
  nrAL=0;
  ncAL=0;
  nrncAL=0;
  nrAU=0;
  ncAU=0;
  nrncAU=0;
  j=0; // true column index
  cntL=cntU=0;
#ifdef PRINT_CHECK
  for (m=0; m<n; m++) {
      if (idxposL[m]) {
	 printf("1.idxposL[%ld]=%ld\n",m,idxposL[m]);
	 fflush(stdout);
      } // end if
  } // end for m	
  for (m=0; m<n; m++) {
      if (idxposU[m]) {
	 printf("1.idxposU[%ld]=%ld\n",m,idxposU[m]);
	 fflush(stdout);
      } // end if
  } // end for m	
#endif
  for (i=0; i<nblocks; i++) {
      // store current block size for easy access
      bi=blocksize[i];

      r=j; s=r+bi;
      for (k=0; k<bi; k++) {
	  // check sub-/super block diagonal entries of block i
	  jj=p[j]; // shortcut
	  for (m=0; m<A->ncol[jj]; m++) {
	      t=invq[A->rowind[jj][m]];
	      // t is below the diagonal block
	      if (t>=s) {
		 // new nonzero row index found in this block column
		 if (!idxposL[t]) {
		    idxL[cntL]=t;
		    idxposL[t]=++cntL;
		    (BL->nblockrow[i])++;
		 } // end if
	      } // end if t>=s      
	      // t is above the diagonal block
	      else if (t<r) {
		 i2=invblock[t];
		 // new nonzero row index found in this block column
		 if (!idxposU[i2]) {
		    idxU[cntU]=i2;
		    idxposU[i2]=++cntU;
		    (BUT->nblockrow[i2])++;
		 } // end if
	      } // end if t<r
	  } // end for m
	  for (m=0; m<cntU; m++) {
	      t=idxU[m];
	      idxposU[t]=0;
	  } // end for m
	  cntU=0;
#ifdef PRINT_CHECK
	  for (m=0; m<n; m++) {
	      if (idxposU[m]) {
		 printf("2.idxposU[%ld]=%ld\n",m,idxposU[m]);
		 fflush(stdout);
	      } // end if
	  } // end for m	
#endif
	  j++;
      } // end for k
      // reset idxposL
      for (m=0; m<cntL; m++) {
	  t=idxL[m];
	  idxposL[t]=0;
      } // end for m
      cntL=0;
#ifdef PRINT_CHECK
      for (m=0; m<n; m++) {
	  if (idxposL[m]) {
	     printf("3.idxposL[%ld]=%ld\n",m,idxposL[m]);
	     fflush(stdout);
	  } // end if
      } // end for m	
#endif
      // maximum number of sub-diagonal indices
      nrAL=MAX(nrAL,BL->nblockrow[i]);
      // maximum block size sub-diagonal block
      nrncAL=MAX(nrncAL,bi*BL->nblockrow[i]);
  } // end for i
  ncAD=0;
  for (i=0; i<nblocks; i++) {
      // printf("%ld,%ld\n",BL->nblockrow[i],BUT->nblockrow[i]);fflush(stdout);
      bi=blocksize[i];
      ncAD=MAX(ncAD,bi);
      // maximum number of super-diagonal indices
      nrAU=MAX(nrAU,BUT->nblockrow[i]);
      // maximum block size super-diagonal block
      nrncAU=MAX(nrncAU,bi*BUT->nblockrow[i]);
  } // end for i
  // printf("%ld,%ld,%ld,%ld\n",nrncAL,nrncAU,nrAL,nrAU);fflush(stdout);
  // finally we have a relatively precise estimate how big the buffers should be
  AvalD     =(FLOAT  *) malloc((size_t)ncAD*ncAD*sizeof(FLOAT));
  AvalL     =(FLOAT *)  malloc((size_t)nrncAL   *sizeof(FLOAT));
  AvalU     =(FLOAT *)  malloc((size_t)nrncAU   *sizeof(FLOAT));
  Arowind   =(integer *)malloc((size_t)nrAL     *sizeof(integer));
  Acolind   =(integer *)malloc((size_t)nrAU     *sizeof(integer));
#ifdef _PROFILING_
  time_init_matrix+=omp_get_wtime()-timeBegin;
#endif

  
  // **************************************************************************
  // *****                       main loop                                *****
  // **************************************************************************
  i=0; // index for the blocks of A
  j=0; // index for the scalar columns/rows of A
  i_bldu=0; // index for the block of the ILU, 0<=i_bldu<=i
  cntL=cntU=0; // counter for temporary index lists
#ifdef PRINT_CHECK
  for (m=0; m<n; m++) {
    if (idxposL[m]) {
      printf("4.idxposL[%ld]=%ld\n",m,idxposL[m]);
      fflush(stdout);
    } // end if
  } // end for m	
  for (m=0; m<n; m++) {
    if (idxposU[m]) {
      printf("4.idxposU[%ld]=%ld\n",m,idxposU[m]);
      fflush(stdout);
    } // end if
  } // end for m	
#endif
  cnt_fill=0; // counter for additional fill-in when merging blocks
  while (j<n) {

        // store current block size for easy access
	bi=blocksize[i];
	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	if (bi>ncAD) {
	   ncAD=bi;
	   AvalD=(FLOAT *)realloc(AvalD,(size_t)ncAD*ncAD*sizeof(FLOAT));
	}
        // clear buffers
        pD=AvalD;
	for (k=0; k<bi*bi; k++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
	    pD[k]=0.0;
#else
	    pD[k].r=pD[k].i=0.0;
#endif
	}
	for (k=0; k<nrncAL; k++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
	    AvalL[k]=0.0;
#else
	    AvalL[k].r=AvalL[k].i=0.0;
#endif
	}
	for (k=0; k<nrncAU; k++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
	    AvalU[k]=0.0;
#else
	    AvalU[k].r=AvalU[k].i=0.0;
#endif
	}
	// for (k=0; k<nrAL; k++) 
	//    Arowind[k]=0;
	// for (k=0; k<nrAU; k++) 
	//    Acolind[k]=0;
#ifdef _PROFILING_
	time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	        

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif    
	// ----------------------------------------------
	// ----------------------------------------------
	// ------- extract COLUMNS of A of block i ------
	// pointer to the current column inside diagonal block i
        pD=AvalD;
	// pointer to the current strict lower triangular block
	pL=AvalL;
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 start
        r=j; s=r+bi;

	// prepare idxL to carry the diagonal block indices in increasing
	// order
	for (m=r; m<s; m++) {
	    idxL[cntL]=m;
	    idxposL[m]=++cntL;
	} // end for m
	
	// column indices are r:s-1
	BL->colind[i_bldu] =(integer *)malloc((size_t)bi*sizeof(integer));
	BiD->colind[i_bldu]=(integer *)malloc((size_t)bi*sizeof(integer));
	BUT->colind[i_bldu]=(integer *)malloc((size_t)bi*sizeof(integer));
	BL->nblockcol[i_bldu]=BiD->nblockcol[i_bldu]=BUT->nblockcol[i_bldu]=bi;
	// number of extracted rows
	nrowind=0;
	// loop through the columns r:s-1 inside block column i
	for (k=0; k<bi; k++,pD+=bi,pL+=BL->nblockrow[i_bldu]) {
	    // trivial column indices 
	    BL->colind[i_bldu][k]=BiD->colind[i_bldu][k]=BUT->colind[i_bldu][k]=j;

	    jj=p[j]; // shortcut
	    // check COLUMN j of A(q,p)
	    ncol=A->ncol[jj];
	    rowind=A->rowind[jj];
	    for (m=Afirst[j]; m<ncol; m++) {
 	        // row index tt of A(tt,p[j])
	        tt=rowind[m];
		// tt=q[t]
	        t=invq[tt];
		
		// extract associated numerical value
		val=A->val[jj][m];
		if (SL!=NULL) {
		   rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		   val*=rval;
#else
		   val.r*=rval;
		   val.i*=rval;
#endif
		}
		if (SR!=NULL) {
		   rval=SR[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		   val*=rval;
#else
		   val.r*=rval;
		   val.i*=rval;
#endif
		}
		
		// t is inside the diagonal block
		if (r<=t && t<s) {
		   pD[t-r]=val;
		} // end if r<=t<s
		// t is below the diagonal block
		else if (t>=s) {
		   // check whether this is a new row or not
		   tt=idxposL[t];
		   // new nonzero row index found in this column
		   if (!tt) {
 		      idxL[cntL]=t;       // bookmark index and ...
		      idxposL[t]=++cntL;  // its physical position, shifted by 1
		      tt=cntL;            // return its assigned location       
		      // if (nrowind>=nrAL) {
		      //    nrAL+=10;
		      //    Arowind=(integer *)realloc(Arowind,(size_t)nrAL*sizeof(integer));
		      // }
		      // bookmark index for later purposes
		      // Arowind[nrowind++]=t;
		      nrowind++;
		      // nrowind and cntL are identical; similarly, idxL and Arowind coincide
		   } // end if
		   // whoops, our estimate about the number of nonzeros in the current block row was
		   // too small (no idea how this could have happened)
		   if (nrowind>BL->nblockrow[i]) {
		      // increase number of off-diagonal row indices in this block
		      (BL->nblockrow[i])++;
		      nrncAL=MAX(nrncAL,bi*BL->nblockrow[i]);
		      // if necessary, realloc memory
		      AvalL=(FLOAT *)realloc(AvalL,(size_t)nrncAL*sizeof(FLOAT));
		      // pointer to the last existing element
		      pL=AvalL+(k+1)*(BL->nblockrow[i]-1)-1;
		      // shift already inserted old entries
		      for (jj=k; jj>0; jj--) {
			  // set additional element in column k to 0
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			  pL[jj+1]=0.0;
#else
			  pL[jj+1].r=pL[jj+1].i=0.0;
#endif
			  // shift old column k
			  for (mm=0; mm<(BL->nblockrow[i]-1); mm++,pL--) {
			      pL[jj]=*pL;
			  }
		      }
		      // set additional element in column 0 to 0
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		      pL[1]=0.0;
#else
		      pL[1].r=pL[1].i=0.0;
#endif
		      // re-adjust pointer w.r.t. new number of rows
		      pL=AvalL+k*BL->nblockrow[i];
		   } // end if
		   // insert sub-diagonal entry
		   tt-=bi; // remember that idxL, idxposL include diagonal block
		   pL[tt-1]=val;
		} // end else-if t>=s      
	    } // end for m
	    j++;
	} // end for k

	/*
	// reset idxposL
	for (m=0; m<cntL; m++) {
	    t=idxL[m];
	    idxposL[t]=0;
	} // end for m
	cntL=0;
	*/
	// ----- end extract COLUMNS of A of block i ----
	// ----------------------------------------------
	// ----------------------------------------------

	

	// ----------------------------------------------
	// ----------------------------------------------
        // -------- extract ROWS of A of block i --------
	// pointer to the current strict upper triangular block
	pU=AvalU;
	// reset j
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 start
        j=r; s=r+bi;
	// number of extracted columns
	ncolind=0;
	// loop through the rows r:s-1 inside block column i and search
	// for the required columns
	for (k=0; k<bi; k++,pU+=BUT->nblockrow[i_bldu]) {
	    jj=p[j]; // shortcut
	    // check ROW j of A(q,p)	    
	    // since A is stored by columns we access row j via the linked list
	    jt=Ahead[j];
	    while (jt>=0) {
	          // printf("scanning column %ld\n",jt);
	          jjt=p[jt]; // shortcut
		  // check COLUMN jt of A(q,p), its leading entry must be j
	          m=Afirst[jt];
		  // row index tt of A(tt,p[jt])
		  tt=A->rowind[jjt][m];
		  // tt=q[t], we must have j=t, i.e. we are considerung A(q[j],p[jt])
		  t=invq[tt];
		  
		  // extract associated numerical value
		  val=A->val[jjt][m];
		  if (SL!=NULL) {
		     rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     val*=rval;
#else
		     val.r*=rval;
		     val.i*=rval;
#endif
		  }
		  if (SR!=NULL) {
		     rval=SR[jjt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		     val*=rval;
#else
		     val.r*=rval;
		     val.i*=rval;
#endif
		  }
		  
		  // jt is right located from the diagonal block
		  if (jt>=s) {
		     // check whether this is a new column or not
		     tt=idxposU[jt];
		     // new nonzero column index found in this block row
		     if (!tt) {
		        idxU[cntU]=jt;      // bookmark index and ...		
			idxposU[jt]=++cntU; // its physical position, shifted by 1
			tt=cntU;	    // return its assigned location       
			// if (ncolind>=nrAU) {
			//   nrAU+=10;
			//   Acolind=(integer *)realloc(Acolind,(size_t)nrAU*sizeof(integer));
			// }
			// bookmark index for later purposes
			// Acolind[ncolind++]=jt;
			ncolind++;
			// ncolind and cntU are identical; similarly, idxU and Acolind coincide
		     } // end if
		     // whoops, our estimate about the number of nonzeros in the current (transposed) block row was too small
		     if (ncolind>BUT->nblockrow[i]) {
		        // increase number of off-diagonal row indices in this block
		        (BUT->nblockrow[i])++;
			nrncAU=MAX(nrncAU,bi*BUT->nblockrow[i]);
			// if necessary, realloc memory
			AvalU=(FLOAT *)realloc(AvalU,(size_t)nrncAU*sizeof(FLOAT));
			// pointer to the last existing element
			pU=AvalU+(k+1)*(BUT->nblockrow[i]-1)-1;
			// shift already inserted old entries
			for (jj=k; jj>0; jj--) {
			    // set additional element in column k to 0
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			    pU[jj+1]=0.0;
#else
			    pU[jj+1].r=pU[jj+1].i=0.0;
#endif
			    // shift old column k
			    for (mm=0; mm<(BUT->nblockrow[i]-1); mm++,pU--) {
			        pU[jj]=*pU;
			    } // end for mm
			} // end for jj
			// set additional element in column 0 to 0
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			pU[1]=0.0;
#else
			pU[1].r=pU[1].i=0.0;
#endif
			// re-adjust pointer w.r.t. new number of rows
			pU=AvalU+k*BUT->nblockrow[i];
		     } // end if
		     // insert right-diagonal entry
		     pU[tt-1]=val;
		  } // end else-if jt>=s      
	          jt=Alist[jt];
	    } // end while jt>=0
	    // ---------------------------------------------------------------
	    // ----- update linked list for the next row  greater than j -----
	    // beginning of the linked list
	    jt=Ahead[j];
	    // while list is non-empty
	    while (jt>=0) {
	          // position of the current leading entry in column jt referring to A(q[j:end],p[jt])
	          m=Afirst[jt];
		  jjt=p[jt]; // shortcut
		  // row index tt of A(tt,p[jt])
		  tt=A->rowind[jjt][m];
		  // q[t]=tt,  we should have t=j
		  t=invq[tt];
		  // increment position of the next leading entry in column jt
		  Afirst[jt]=++m;
		  // buffer next column to be considered
		  jt_next=Alist[jt];
		  // we have to make sure that we are still inside column jt
		  if (m<A->ncol[jjt]) {
		     // row index tt of A(tt,p[jt]) we should have t=j
		     tt=A->rowind[jjt][m];
		     // q[t]=tt
		     t=invq[tt];
		     // add new entry to the head of the list
		     Alist[jt]=Ahead[t];
		     Ahead[t]=jt;
		  } // end if 
		  // restore next column
		  jt=jt_next;
	    } // end while jt>=0      
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
	    
	    // advance row counter
	    j++;
	} // end for k

	/*
	// reset idxposU
	for (m=0; m<cntU; m++) {
	    t=idxU[m];
	    idxposU[t]=0;
	} // end for m
	cntU=0;
	*/
        // ------ end extract ROWS of A of block i ------
	// ----------------------------------------------
	// ----------------------------------------------

	/*
	if (nrowind!=cntL-bi) {
	   printf("step %ld, nrowind=%ld, cntL=%ld, bi=%ld\n",
		  i_bldu,nrowind,cntL,bi);
	   fflush(stdout);
	}
	if (ncolind!=cntU) {
	   printf("step %ld, ncolind=%ld, cntU=%ld\n",
		  i_bldu,ncolind,cntU);
	   fflush(stdout);
	}
	*/
	// our estimate for BL->nblockrow[i] was too big, no idea
	// how this could have happened, maybe due to progressive
	// aggregation
	if (BL->nblockrow[i]>nrowind) {
           pL=AvalL+nrowind;
	   pD=AvalL+BL->nblockrow[i];
           for (k=1; k<bi; k++,pD+=BL->nblockrow[i]-nrowind) {
	       for (m=0; m<nrowind; m++) {
		   *pL++=*pD++;
	       } // end for m
           } // end for k
	   BL->nblockrow[i]=nrowind;
        } // end if
	// our estimate for BUT->nblockrow[i] was too big, no idea
	// how this could have happened, maybe due to progressive
	// aggregation
	if (BUT->nblockrow[i]>ncolind) {
           pU=AvalU+ncolind;
	   pD=AvalU+BUT->nblockrow[i];
           for (k=1; k<bi; k++,pD+=BUT->nblockrow[i]-ncolind) {
	       for (m=0; m<ncolind; m++) {
		   *pU++=*pD++;
	       } // end for m
           } // end for k
	   BUT->nblockrow[i]=ncolind;
        } // end if
	

	// sort AvalL and AvalU such that their row indices are sorted
	// in increasing order
	// QSORTGNL(AvalL,Arowind,idxL,&nrowind,&bi);
	QSORTGNL(AvalL,idxL+bi,idxbuff,&nrowind,&bi);
	// re-init off-diagonal positions, now ordered increasingly
	for (m=bi; m<cntL; ) {
	    t=idxL[m];
	    idxposL[t]=++m;
	} // end for m
	// copy off-diagonal idxL indices to Arowind
	if (nrowind>=nrAL) {
	   nrAL=MAX(nrowind,nrAL+10);
	   Arowind=(integer *)realloc(Arowind,(size_t)nrAL*sizeof(integer));
	} // end if
	memcpy(Arowind,idxL+bi,nrowind*sizeof(integer));
	
	// QSORTGNL(AvalU,Acolind,idxU,&ncolind,&bi);
	QSORTGNL(AvalU,idxU,idxbuff,&ncolind,&bi);
	// copy off-diagonal idxL indices to Arowind
	if (ncolind>=nrAU) {
	   nrAU=MAX(ncolind,nrAU+10);
	   Acolind=(integer *)realloc(Acolind,(size_t)nrAU*sizeof(integer));
	} // end if
	memcpy(Acolind,idxU,ncolind*sizeof(integer));
	

	// for testing, print block structure as extracted from A
#ifdef PRINT_INFO0
	printf("block %3ld of A\n",i);
	printf("diagonal block of A\n        ");
	for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	    printf("%8ld",BiD->colind[i_bldu][m]);
	printf("\n");
	for (k=0; k<BiD->nblockcol[i_bldu]; k++) {
	    printf("%8ld",BiD->colind[i_bldu][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",AvalD[m*BiD->nblockcol[i_bldu]+k]);
	    printf("\n");
#else	    
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",AvalD[m*BiD->nblockcol[i_bldu]+k].r);
	    printf("\n        ");
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",AvalD[m*BiD->nblockcol[i_bldu]+k].i);
	    printf("\n");
#endif
	}
	printf("\n");

	printf("sub-diagonal block of A\n        ");
	if (nrowind>0) {
	   for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	       printf("%8ld",BiD->colind[i_bldu][m]);
	   printf("\n");
	   for (k=0; k<nrowind; k++) {
	       printf("%8ld",idxL[bi+k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	           printf("%8.1le",AvalL[m*nrowind+k]);
#else
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",AvalL[m*nrowind+k].r);
	       printf("\n        ");
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",AvalL[m*nrowind+k].i);
#endif
	       printf("\n");
	   } // end for k
	}
	printf("\n");

	printf("transposed super-diagonal block of A\n        ");
	if (ncolind>0) {
	   for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	       printf("%8ld",BiD->colind[i_bldu][m]);
	   printf("\n");
	   for (k=0; k<ncolind; k++) {
	       printf("%8ld",idxU[k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",AvalU[m*ncolind+k]);
#else
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",AvalU[m*ncolind+k].r);
	       printf("\n        ");
	       for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",AvalU[m*ncolind+k].i);
#endif
	       printf("\n");
	   } // end for k
	}
	printf("\n");
#endif
	// reset idxposU
	for (m=0; m<cntU; m++) {
	    t=idxU[m];
	    idxposU[t]=0;
	} // end for m
	cntU=0;
#ifdef PRINT_CHECK
	for (m=0; m<n; m++) {
	    if (idxposU[m]) {
	       printf("5.idxposU[%ld]=%ld\n",m,idxposU[m]);
	       fflush(stdout);
	    } // end if
	} // end for m	
#endif
	
#ifdef _PROFILING_
	time_extract_block+=omp_get_wtime()-timeBegin;
#endif



	// now that block column/row i has been extracted into three blocks
	// we can compute block column i_bldu of L as well as block row
	// i_bldu of U
	// To do so, we use linked list for rows of BL associated with the
	// block columns of BL (similarly for BUT)
	

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// ----------------------------------------------------------------------
	// ----- compute block column i_bldu of BL using linked list of BUT -----
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 of A start
        r=BiD->colind[i_bldu][0]; s=r+bi;


	// first pass: set up nonzero row structure of D+L before dropping
	
	// beginning of the linked list for block column i
	// this information is obtained from the sparsity pattern of U
	it=BUThead[i_bldu];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning block column %ld of BL\n",it); fflush(stdout);
#endif
	      // check block COLUMN it of BL
	      rowind=BL->rowind[it];
	      nrow=BL->nblockrow[it];
	      if (cntL>bi) {
		 // idxL includes more indices than just the diagonal block,
		 for (m=BLfirst[it]; m<nrow; m++) {
		     t=rowind[m];
		     // skip entries above or in the diagonal block
		     if (t>=s)
		        break;
		 } // end for m
		 // use idxU temporarily for the new indices
		 for (; m<nrow; m++) {
		     t=rowind[m];
		     // only store new indices
		     if (!idxposL[t]) {
		        idxU[cntU++]=t;
			// printf("new fill (row %ld,bcol %ld->bcol %ld)\n", t,it,i_bldu);fflush(stdout);
		     } // new fill
		 } // end for m
		 // now idxL[0,...,cntL-1] and idxU[0,...,cntU-1] are
		 // distinct index arrays each with indices in increasing order
		 // we merge them in idxL
		 jj=cntL-1; ii=cntU-1;
		 i2=cntL+cntU;
		 tt=idxL[jj];
		 t=idxU[ii]; 
		 while (ii>=0) {
		       // new index greater than given index
		       if (t>tt) {
			  // printf("insert new %ld at the end\n", t);fflush(stdout);
			  idxposL[t]=i2--;
			  idxL[i2]=t;
			  ii--;
			  t=idxU[ii]; 
		       }
		       else { // tt>t
			  // printf("insert old %ld at the end\n", tt);fflush(stdout);
			  while (tt>t) {
			        idxposL[tt]=i2--;
				idxL[i2]=tt;
				jj--;
				tt=idxL[jj];
			  } // end while
		       } // end if-else
		 } // end while
		 cntL+=cntU;
		 cntU=0;
	      }
	      else { // idxL only includes the diagonal block, we can append the indices directly
		 // printf("cntL=%ld,cntU=%ld\n",cntL,cntU);fflush(stdout);
		 for (m=BLfirst[it]; m<nrow; m++) {
		     t=rowind[m];
		     // skip entries above or in the diagonal block
		     if (t>=s) 
		        break;
		 } // end for m
		 for (; m<nrow; m++) {
		     t=rowind[m];
		     // append new indices
		     idxL[cntL]=t;
		     idxposL[t]=++cntL;
		     // printf("new fill (row %ld,bcol %ld->bcol %ld)\n", t,it,i_bldu);fflush(stdout);
		 } // end for m
	      } // end if-else cntL
	      
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		  if (idxposU[m]) {
		     printf("6.idxposU[%ld]=%ld\n",m,idxposU[m]);
		     fflush(stdout);
		  } // end if
	      } // end for m	
#endif

	      // advance to next required block column of BL as requested by BUT
	      it=BUTlist[it];
	} // end while
	// end first pass
#ifdef PRINT_INFO0
	printf("total nonzero row indices column %ld of BL\n",i); fflush(stdout);
	for (m=0; m<cntL; m++) {
  	    printf("%8ld",idxL[m]);
	} // end for m
	printf("\n");fflush(stdout);
#endif

	
	// preparing second pass
	// now we know the size of the current block column before dropping
	// buff is going to store block column i of L+D temporarily
	if (cntL*bi>nbuff) {
	   // increase buffer
	   nbuff=ELBOW_BUFF*(cntL*bi)+BUFF_EXT;
	   buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	} // end if
	// copy diagonal block of A to buff
	pb=buff; pD=AvalD;
	for (m=0; m<bi; m++,pb+=cntL,pD+=bi) {
	    // copy columns r,...,s-1 of the diagonal block
	    memcpy(pb,pD,(size_t)bi*sizeof(FLOAT));
	} // end for m
	// init remaining buff with 0 and scatter sub-diagonal block of A into buff
 	pb=buff; pL=AvalL;
	for (m=0; m<bi; m++,pb+=cntL,pL+=nrowind) {
	    // init rows with 0
	    for (jj=bi; jj<cntL; jj++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	        pb[jj]=0.0;
#else
		pb[jj].r=pb[jj].i=0.0;
#endif
	    } // end for k
	    // copy rows of the sub-diagonal block of A
	    for (jj=0; jj<nrowind; jj++) {
	        // row index t
	        t=Arowind[jj];
		// associated position shifted by 1
		tt=idxposL[t];
		pb[tt-1]=pL[jj];
	    } // end for jj
	} // end for m
#ifdef _PROFILING_
	time_LU_update_pass_1+=omp_get_wtime()-timeBegin;
#endif
	
	// second pass: compute block column i of BL using level-3-BLAS
	// beginning of the linked list of BUT, the nonzero pattern of BUT
	// reflects the requested columns of BL
	it=BUThead[i_bldu];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning again block column %ld of BL\n",it); fflush(stdout);
#endif
	      // level-3-BLAS update C+=AB^T

#ifdef _PROFILING_
	      timeBegin = omp_get_wtime();
#endif
	      // step 1: rows of block column "it" of BUT refer to rows of B
	      //         These can be left in-place for the level-3-BLAS update
	      // read the associated rows from BUT for B
	      jj=BUTfirst[it];        // start of row r and higher
	      ldb=BUT->nblockrow[it]; // leading dimension (total number of rows)
	      flag=-1;
	      rowind=BUT->rowind[it]; // shortcut
	      while (jj<ldb && flag) {
		    // row indices between r and s-1
		    t=rowind[jj];
		    // row indices are inside the current block column i
		    if (t<s) {
		       idxU[cntU]=t;      // bookmark index and ...
		       idxposU[t]=++cntU; // its position
		       jj++;
		    }
		    else // block column i exceeded
		       flag=0;
	      } // end while
	      // now idxU, idxposU refer to the columns of block column i
	      // to be updated with "it"
#ifdef PRINT_INFO0
	      printf("column indices required to update with column %ld of BL\n",it); fflush(stdout);
	      for (m=0; m<cntU; m++) {
		  printf("%8ld",idxU[m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef PRINT_INFO0
	      printf("row indices required to update with column %ld of BL\n",it); fflush(stdout);
	      for (m=BLfirst[it]; m<BL->nblockrow[it]; m++) {
		  printf("%8ld",BL->rowind[it][m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	      time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif

	      l=BL->nblockcol[it];
	      ii=BLfirst[it];        // start of row r and higher in BL
	      lda=BL->nblockrow[it]; // leading dimension of block column "it" of BL
	      ldc=lda-ii;            // number of rows of the block to be updated
	      // if there is anything to update 
	      if (ldc && cntU && l) {
#ifdef _PROFILING_
		 timeBegin=omp_get_wtime();
#endif
		 // single-column-update, do this by hand and skip level-3-BLAS
		 // compress index arrays into one array
		 pi=BL->rowind[it]+ii; // row indices to be considered
		 for (k=0; k<ldc; k++) {
		     // which row indices are required?
		     tt=*pi++;
		     idxbuff[k]=idxposL[tt];
		 } // end for k
#ifdef _PROFILING_
		 time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif
		 if (cntU==1 && l==1) {
#ifdef _PROFILING_
		    timeBegin=omp_get_wtime();
#endif
		    // extract column index t
		    t=idxU[0];
		    // scatter to column t of buff
		    pb=buff+cntL*(t-r)-1;   // associated column of buff,shifted
		    pL=BL->valE[it]  +ii; // numerical values of column BL{it}
		    pU=BUT->valE[it] +BUTfirst[it]; // numerical value BUT{it}
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
		    val=-(*pU);
#else
		    val.r=-pU->r;
		    val.i=-pU->i;
#endif

#ifdef _USE_MKL_
		    // use sparse AXPYI pb[idxbuff]= val * pL + pb[idxbuff]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    CBLAS_AXPYI(ldc, val,  pL,idxbuff, pb); 
#else
		    CBLAS_AXPYI(ldc, &val, pL,idxbuff, pb); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
		    
#else // hand-coded loop
		    for (k=0; k<ldc; k++) {
			// position of i2 in buff
			i2=idxbuff[k];
			// downdate buff - BL{i}*BUT{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			pb[i2]+=(*pL++)*val;
#else
			pb[i2].r+=pL->r*val.r-pL->i*val.i;
			pb[i2].i+=pL->r*val.i+pL->i*val.r;
			pL++;
#endif
		    } // end for k
#endif //-else _USE_MKL_
#ifdef _PROFILING_
		    time_scalar_update+=omp_get_wtime()-timeBegin;
#endif
		    
		 }
		 else if ((double)ldc*cntU<=SMALL_BLOCK &&
			  (double)ldc*l   <=SMALL_BLOCK &&
			  (double)l*cntU  <=SMALL_BLOCK) {
 		    // use hand-coded loops and update in-place

		    // rank-1-update
		    pU=BUT->valE[it] +BUTfirst[it]; // numerical value BUT{it}
		    if (l==1) {
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
 		       for (m=0; m<cntU; m++,pU++) {
			   // numerical values of column BL{it}
			   pL=BL->valE[it]  +ii; 
			   // extract column index t
			   t=idxU[m];
			   // directly scatter to column t of buff
			   pb=buff+cntL*(t-r)-1;   // associated column of buff, shifted
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   val=-(*pU);
#else
			   val.r=-pU->r;
			   val.i=-pU->i;
#endif
			   
			   
#ifdef _USE_MKL_
			   // use sparse AXPYI pb[idxbuff]= val * pL + pb[idxbuff]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   CBLAS_AXPYI(ldc, val,  pL,idxbuff, pb); 
#else
			   CBLAS_AXPYI(ldc, &val, pL,idxbuff, pb); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
			   
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       // downdate buff - BL{i}*BUT{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			       pb[i2]+=(*pL++)*val;
#else
			       pb[i2].r+=pL->r*val.r-pL->i*val.i;
			       pb[i2].i+=pL->r*val.i+pL->i*val.r;
			       pL++;
#endif
			   } // end for k
#endif //-else _USE_MKL_

		       } // end for m
#ifdef _PROFILING_
		       time_rank_1_update+=omp_get_wtime()-timeBegin;
#endif
		    }
		    else { // higher rank using DOT product
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
		       for (m=0; m<cntU; m++,pU++) {
			   // numerical values of column BL{it}
			   pL=BL->valE[it]  +ii; 
			   // extract column index t
			   t=idxU[m];
			   // directly scatter to column t of buff
			   pb=buff+cntL*(t-r)-1;   // associated column of buff, shifted
			   for (k=0; k<ldc; k++) {
			       // position i2 buff
			       i2=idxbuff[k];
			       // downdate buff - BL{i}*BUT{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			       pb[i2]-=DISTDOT(&l, pL,&lda, pU,&ldb);
#else
#ifdef _USE_MKL_			       
			       DISTDOTU(&val, &l, pL,&lda, pU,&ldb);
#else
			       val=DISTDOTU(&l, pL,&lda, pU,&ldb);
#endif
			       pb[i2].r-=val.r;
			       pb[i2].i-=val.i;
#endif
			       pL++;
			   } // end for k
		       } // end for m
#ifdef _PROFILING_
		       time_small_block_update+=omp_get_wtime()-timeBegin;
#endif
		    } // end if-else l=1
		 }
		 else {
#ifdef _PROFILING_
		    timeBegin=omp_get_wtime();
#endif
		    // check whether we are able to bypass the gather/scatter part
		    // this is the case, whenever the selected column indices are
		    // a contiguous sequence and the associated row indices are
		    // contiguous sequence as well
		    if (idxbuff[ldc-1]-idxbuff[0]+1==ldc &&
			idxU[cntU-1]  -idxU[0]   +1==cntU) {
		       skip_gthr_sctr=-1;
		       // printf("gather/scatter could be skipped\n"); fflush(stdout);
		    }
		    else
		       skip_gthr_sctr=0;

 		    // step 2: temporarily gather the data from buff to (buff)C to
		    //         accomodate the level-3-BLAS update
		    if (!skip_gthr_sctr) {
		       // extract the rows of buff associated with the rows of BL
		       if (cntU*ldc>nbuffC) {
			  nbuffC=ELBOW_BUFF*(cntU*ldc)+BUFF_EXT;
			  buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
		       } // end if
	      
		       // gather the associated rows and columns from buff to (buff)C
		       pC=buffC;
		       for (m=0; m<cntU; m++) {
			   // extract column index t
			   t=idxU[m];
			   // gather column t from buff into buffC
			   pb=buff+cntL*(t-r)-1; // associated column of buff, shifted

#ifdef _USE_MKL_
			   CBLAS_GTHR(ldc, pb, pC,idxbuff);
			   pC+=ldc;
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 buff
			       i2=idxbuff[k];
			       *pC++=pb[i2];
			   } // end for k
#endif //-else _USE_MKL_
			
		       } // end for m
		    } // end if NOT skip_gthr_sctr
#ifdef _PROFILING_
		    time_level_3_blas_gather+=omp_get_wtime()-timeBegin;
#endif
		    
		    // now perform level-3-BLAS update
		    // C = -1*A*B^T+1*C
		    transa="n"; transb="t";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    alpha=-1.0;
		    beta = 1.0;
#else
		    alpha.r=-1.0; alpha.i=0.0;
		    beta.r = 1.0; beta.i =0.0;
#endif
		    // level-3 BLAS performed in-place
		    if (skip_gthr_sctr)
		       GEMM(transa,transb, &ldc,&cntU,&l, &alpha,
			    BL->valE[it]+BLfirst[it],&lda,
			    BUT->valE[it]+BUTfirst[it],&ldb, &beta,
			    buff + idxbuff[0]-1 + cntL*(idxU[0]-r),&cntL, 1,1);
		    else // level-3 BLAS uses a cached copy
		       GEMM(transa,transb, &ldc,&cntU,&l, &alpha,
			    BL->valE[it]+BLfirst[it],&lda,
			    BUT->valE[it]+BUTfirst[it],&ldb, &beta,
			    buffC,&ldc, 1,1);
	      
#ifdef _PROFILING_
		    timeBegin2=omp_get_wtime();
#endif
		    if (!skip_gthr_sctr) {
		       // scatter the result back to the associated rows
		       // and columns of buff
		       pC=buffC;
		       for (m=0; m<cntU; m++) {
		           // extract column index t
			   t=idxU[m];
			   // scatter to column t of buff
			   pb=buff+cntL*(t-r)-1; // associated column of buff, shifted

#ifdef _USE_MKL_
			   CBLAS_SCTR(ldc, pC, idxbuff, pb);
			   pC+=ldc;
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       pb[i2]=*pC++;
			   } // end for k
#endif //-else _USE_MKL_
			
		       } // end for m
		    } // end if NOT skip_gthr_sctr
#ifdef _PROFILING_
		    time_level_3_blas_update+=omp_get_wtime()-timeBegin;
		    time_level_3_blas_scatter+=omp_get_wtime()-timeBegin2;
#endif
		 } // end if-elseif-else cntU=1 and l=1
	      } // end if ldc> and cntU>0 and l>0
		 
#ifdef _PROFILING_
	      timeBegin=omp_get_wtime();
#endif
	      // update BUTfirst, this automatically excludes the diagonal
	      // block when block ROW i is updated
	      BUTfirst[it]=jj;

	      // reset idxposU
	      for (m=0; m<cntU; m++) {
		  t=idxU[m];
		  idxposU[t]=0;
	      } // end for m
	      cntU=0;
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		  if (idxposU[m]) {
		     printf("7.idxposU[%ld]=%ld\n",m,idxposU[m]);
		     fflush(stdout);
		  } // end if
	      } // end for m	
#endif
	      
	      // advance to the next block column
	      it=BUTlist[it];
#ifdef _PROFILING_
	      time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	} // end while
	// end second pass

	// scalar column, skip LAPACK
	if (bi==1) {
	   if (!invert_blocks)
	      pivots[0]=1;
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // the diagonal entry is buff[0]
  	   rval=FABS2(buff[0]); // |diagonal_entry|^2
	   
	   // zero pivot
	   flag=0;
	   as=rs=0.0;
	   if (rval==0.0) {
	      ierr=1;
	      flag=1;
	      // small absolute shift for a zero pivot
	      as=ABS_THRESHOLD*Amax;
	   }
	   else
	      ierr=0;

	   // tiny pivot?
	   if (!flag && perturbation) {
	      // almost zero pivot
	      if (rval<=macheps*Amax) {
		 flag=1;
		 // small absolute shift for an almost zero pivot
		 as=ABS_THRESHOLD*Amax;
	      }
	      // compute largest sub-diagonal entry in modulus
	      ldc=cntL-1;
	      l=1;
	      ii=I_AMAX(&ldc,buff+1,&l)-1;
	      Ajmax=FABS(buff[ii+1]);
	      if (ABS_THRESHOLD*Ajmax>=rval) {
		 flag=1;
		 // small relative shift for a tiny pivot
		 rs=droptol*REL_THRESHOLD*rval;
	      }
	   } // end if flag and perturbation
	   
	   // lift zero/tiny pivot
	   if (flag && perturbation) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      if (buff[0]>=0) 
		 buff[0]+=as+rs;
	      else
		 buff[0]-=as+rs;
#else
	      if (buff[0].r>=0) 
		 buff[0].r+=as+rs;
	      else
		 buff[0].r-=as+rs;
#endif
	      rval=FABS2(buff[0]); // recompute |diagonal_entry|^2
	      ierr=0;
	   } // end if flag and perturbation

	   if (ierr) {
	      printf("BILU encountered zero pivot\n");

	      // sort original matrix back such that row indices and their
	      // values are taken in increasing order
	      for (j=0; j<n; j++)
		  QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
	      
	      idxL--;idxL=FREE(idxL);
	      idxU--;idxU=FREE(idxU);
	      idxbuff =FREE(idxbuff);
	      idxposL =FREE(idxposL);
	      idxposU =FREE(idxposU);
	      invblock=FREE(invblock);
	      Ahead   =FREE(Ahead);
	      Alist   =FREE(Alist);
	      Afirst  =FREE(Afirst);
	      BLhead  =FREE(BLhead);
	      BLlist  =FREE(BLlist);
	      BLfirst =FREE(BLfirst);
	      BUThead =FREE(BUThead);
	      BUTlist =FREE(BUTlist);
	      BUTfirst=FREE(BUTfirst);
	      AvalD   =FREE(AvalD);
	      AvalL   =FREE(AvalL);
	      AvalU   =FREE(AvalU);
	      Arowind =FREE(Arowind);
	      Acolind =FREE(Acolind);
	      buffC   =FREE(buffC);
	      buff    =FREE(buff);
	      
	      for (i=0; i<nblocks; i++) {
		  BL->colind[i] =FREE(BL->colind[i] );
		  BL->rowind[i] =FREE(BL->rowind[i] );
		  BUT->colind[i]=FREE(BUT->colind[i]);
		  BUT->rowind[i]=FREE(BUT->rowind[i]);
		  BiD->colind[i]=FREE(BiD->colind[i]);
		  BL->valE[i]   =FREE(BL->valE[i]   );
		  BUT->valE[i]  =FREE(BUT->valE[i]  );
		  BiD->valD[i]  =FREE(BiD->valD[i]  );
	      }
	      BL->nblockcol =FREE(BL->nblockcol);
	      BL->nblockrow =FREE(BL->nblockrow);
	      BL->colind    =FREE(BL->colind);
	      BL->rowind    =FREE(BL->rowind);
	      BL->valE      =FREE(BL->valE);
	      
	      BUT->nblockcol=FREE(BUT->nblockcol);
	      BUT->nblockrow=FREE(BUT->nblockrow);
	      BUT->colind   =FREE(BUT->colind   );
	      BUT->rowind   =FREE(BUT->rowind   );
	      BUT->valE     =FREE(BUT->valE     );
	      
	      BiD->nblockcol=FREE(BiD->nblockcol);
	      BiD->colind   =FREE(BiD->colind   );
	      BiD->valD     =FREE(BiD->valD     );
	      
	      BL->nc=BL->nr=BUT->nc=BUT->nr=BiD->nc=BiD->nr=0;
	      BL->nblocks=BUT->nblocks=BiD->nblocks=0;
	      
	      if (flag_ilu1t) {
		 blocksize=FREE(blocksize);
		 blocksize=myblocksize;
		 nblocks  =mynblocks;
	      } // end if flag_ilu1t

	      if (perturbation)
		 buffDGLBLCK=FREE(buffDGLBLCK);
	      
	      return (ierr);
	   } // end if ierr

	   
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	   // logarithm of a real number
	   rval=buff[0];
	   if (rval>=0.0) {
	      locdet.r=log(rval);
	      locdet.i=0.0;
	   }	      
	   else { // log of negative numbers is complex (main branch)
	      locdet.r=log(-rval); 
	      locdet.i=M_PI;
	   }
#else
	   // complex logarithm (main branch)
	   rval=FABS(buff[0]);
	   // 1. real part
	   locdet.r=log(rval);
	   // 1. imaginary part
	   // compute argument of buff[0]
	   if (buff[0].i>=0) 
	      locdet.i=acos(buff[0].r/rval);
	   else
	      locdet.i=-acos(buff[0].r/rval);
#endif

	   
	   rval=FABS2(buff[0]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	   buff[0]=1.0/buff[0];
#else
	   buff[0].r/= rval;
	   buff[0].i/=-rval;
#endif

	   // extract inverse diagonal block to BiD
	   BiD->valD[i_bldu]=(FLOAT *)malloc(sizeof(FLOAT));
	   BiD->valD[i_bldu][0]=buff[0];
#ifdef _PROFILING_
	   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // apply dropping to the strict lower triangular block of BL
	   ldc=cntL-1; // number of rows before dropping
	   rval=sqrt(rval)*droptol; // |diagonal_entry|*droptol
	   pb=buff+1; // start strict lower triangular part
	   l=0;
	   for (ii=0; ii<ldc; ii++,pb++)
	       if (FABS(*pb)>=rval)
		  l++;
	   // now we know the fill
	   BL->nblockrow[i_bldu]=l;
	   BL->valE[i_bldu]=(FLOAT *)malloc((size_t)l*sizeof(FLOAT));
	   BL->rowind[i_bldu]=(integer *)malloc((size_t)l*sizeof(integer));
	   	   
	   // multiply strict lower triangular part by the inverse diagonal entry
	   // and copy it to BL
	   pb=buff+1;
	   pL=BL->valE[i_bldu];
	   pi=BL->rowind[i_bldu];
	   for (ii=0; ii<ldc; ii++,pb++)
	       if (FABS(*pb)>=rval) {
		  *pi++=idxL[ii+1]; // shift idxL by 1 because of the diagonal index
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		  
		  *pL++=(*pb)*(*buff);
#else
		  pL->r=pb->r*buff->r-pb->i*buff->i;
		  pL->i=pb->r*buff->i+pb->i*buff->r;
		  pL++;
#endif		  
	       }
#ifdef _PROFILING_
	       time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	}
	else { // bi>1, block case
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // check prior to the factorizatiion whether there is a column or
	   // row in the diagonal block with a tiny norm
	   if (perturbation) {
	          
	      // check columns of the diagonal block
	      l=1;
	      pb=buff;
	      // small absolute shift for an almost zero row
	      as=ABS_THRESHOLD*Amax;
	      for (jj=0; jj<bi; jj++,pb+=cntL) {
		  // largest entry in column jj, diagonal block
		  ii=I_AMAX(&bi,pb,&l)-1;
		  val=pb[ii];
		  rval=FABS(val);
		  // almost zero column
		  if (rval<=macheps*Amax) {
		     // check if we could shift the diagonal entry instead
		     if (ii!=jj) {
		        v =pb[jj];
		        rv=FABS(v);
		        if (rv>=0.5*rval) {
			   ii=jj;
			   rval=rv;
			   val =v;
			}
		     } // end if ii!=jj
		     // small absolute shift for an almost zero column
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     if (val>=0) 
		        val+=as;
		     else
		        val-=as;
#else
		     if (val.r>=0) 
		        val.r+=as;
		     else
		        val.r-=as;
#endif
		     pb[ii]=val;
		  } // end if
	      } // end for jj
	      // check rows of the diagonal block
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb++) {
		  // largest entry in row jj, diagonal block
		  ii=I_AMAX(&bi,pb,&cntL)-1;
		  val=pb[ii*cntL];
		  rval=FABS(val);
		  // almost zero row
		  if (rval<=macheps*Amax) {
		     // check if we could shift the diagonal entry instead
		     if (ii!=jj) {
		        v =pb[jj*cntL];
		        rv=FABS(v);
		        if (rv>=0.5*rval) {
			   ii=jj;
			   rval=rv;
			   val =v;
			} 
		     } // end if ii!=jj
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     if (val>=0) 
		        val+=as;
		     else
		        val-=as;
#else
		     if (val.r>=0) 
		        val.r+=as;
		     else
		        val.r-=as;
#endif
		     pb[ii*cntL]=val;
		  } // end if
	      } // end for jj
	      // use buffC as work array
	      l=bi*bi;
	      if (l>nbuffDGLBLCK) {
		 nbuffDGLBLCK=ELBOW_BUFF*l+BUFF_EXT;
		 buffDGLBLCK=(FLOAT *)realloc(buffDGLBLCK,
					      (size_t)nbuffDGLBLCK*sizeof(FLOAT));
	      } // end if
	      // make a copy of the diagonal block
	      pb=buff;
	      pC=buffDGLBLCK;
	      ii=1;
	      Anrm=0.0;
	      for (jj=0; jj<bi; jj++,pb+=cntL,pC+=bi) {
		  rv=ASUM(&bi, pb,&ii);
		  Anrm=MAX(Anrm,rv);
		  COPY(&bi, pb,&ii, pC,&ii);
	      } // end for jj
	   } // end if perturbation

	   
	   // factorize and invert diagonal block using LAPACK
	   // use idxU as buffer for pivot indices
	   GETRF(&bi,&bi, buff,&cntL, idxU, &ierr);

	   // check condition number
	   if (ierr || perturbation) {
	      // flag for singular or ill-conditioned system
	      flag=0;
	      // singular case
	      if (ierr) {
		 flag=1;
		 // printf("singular system\n");fflush(stdout);
	      }
	      
	      // if no error occured (i.e. system is non-singular) use LU
	      // decomposition to estimate the condition number of the system
	      if (!flag) {
		 // use buffC as work array
		 l=4*bi;
		 if (l>nbuffC) {
		    nbuffC=ELBOW_BUFF*l+BUFF_EXT;
		    buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
		 } // end if
		 // printf("estimate condition number\n");fflush(stdout);
#if defined _SINGLE_REAL_|| defined _DOUBLE_REAL_
		 GECON("1", &bi, buff, &cntL, &Anrm, &rv, buffC, idxbuff,
		       &jj,1);
#else
		 GECON("1", &bi, buff, &cntL, &Anrm, &rv, buffC+2*bi,
		       (REALS *)buffC, &jj,1);
#endif
		 if (jj) {
		    printf("the %ld-th argument had an illegal value\n",-jj);
		    fflush(stdout);
		    exit (0);
		 }
		 // check whether the system ill-conditioned
		 if (rv<=macheps) {
		    flag=1;
		    // printf("large condition number=1/%8.1le\n",rv);fflush(stdout);
		 }
	      } // end if !flag
	   } // end if ierr or perturbation

	   // if the diagonal block is ill-conditioned or singular
	   if (flag && perturbation) {
	      // restore original entries
	      pb=buff;
	      pC=buffDGLBLCK;
	      ii=1;
	      for (jj=0; jj<bi; jj++,pb+=cntL,pC+=bi)
		  COPY(&bi, pC,&ii, pb,&ii);

	      /*
	      printf("restored system\n"); fflush(stdout);
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb++) {
		  for (ii=0; ii<bi; ii++) {
		      printf("%12.4le",pb[ii*cntL]);
		  }
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
	      */
	      
	      // perturb largest entries in each column slightly
	      as=ABS_THRESHOLD*Amax;
	      l=1;
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb+=cntL) {
		  // largest entry in column jj, diagonal block
		  ii=I_AMAX(&bi,pb,&l)-1;
		  val=pb[ii];
		  rval=FABS(val);
		  rs=droptol*REL_THRESHOLD*rval;
		  // check if we could shift the diagonal entry instead
		  if (ii!=jj) {
		     v =pb[jj];
		     rv=FABS(v);
		     if (rv>=0.5*rval) {
		        ii=jj;
			rval=rv;
			val =v;
			rs=droptol*REL_THRESHOLD*rval;
		     }
		  } // end if ii!=jj
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  if (val>=0) 
		     val+=as+rs;
		  else
		     val-=as+rs;
#else
		  if (val.r>=0) 
		     val.r+=as+rs;
		  else
		     val.r-=as+rs;
#endif
		  pb[ii]=val;
	      } // end for jj
	      // perturb largest entries in each row slightly
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb++) {
		  // largest entry in column jj, diagonal block
		  ii=I_AMAX(&bi,pb,&cntL)-1;
		  val=pb[ii*cntL];
		  rval=FABS(val);
		  rs=droptol*REL_THRESHOLD*rval;
		  // check if we could shift the diagonal entry instead
		  if (ii!=jj) {
		     v =pb[jj*cntL];
		     rv=FABS(v);
		     if (rv>=0.5*rval) {
		        ii=jj;
			rval=rv;
			val =v;
			rs=droptol*REL_THRESHOLD*rval;
		     }
		  } // end if ii!=jj
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  if (val>=0) 
		     val+=as+rs;
		  else
		     val-=as+rs;
#else
		  if (val.r>=0) 
		     val.r+=as+rs;
		  else
		     val.r-=as+rs;
#endif
		  pb[ii*cntL]=val;
	      } // end for jj
	      // for safety, also perturb diagonal entries slightly
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb+=cntL) {
		  val=pb[jj];
		  rval=FABS(val);
		  rs=droptol*REL_THRESHOLD*rval;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  if (val>=0) 
		     val+=as+rs;
		  else
		     val-=as+rs;
#else
		  if (val.r>=0) 
		     val.r+=as+rs;
		  else
		     val.r-=as+rs;
#endif
		  pb[jj]=val;
	      } // end for jj

	      
	      /*
	      printf("shifted system\n"); fflush(stdout);
	      pb=buff;
	      for (jj=0; jj<bi; jj++,pb++) {
		  for (ii=0; ii<bi; ii++) {
		      printf("%12.4le",pb[ii*cntL]);
		  }
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
	      */
	      
	      // recompute LU decomposition
	      GETRF(&bi,&bi, buff,&cntL, idxU, &ierr);
	   } // end if flag and perturbation


	   
	   if (ierr) {
	      if (ierr<0)
		 printf("LAPACK's GETRF: %ld-th argument had an illegal value\n", -ierr);
	      else
		 printf("LAPACK's GETRF routine encountered zero column in step %ld\n", ierr);

	      // sort original matrix back such that row indices and their
	      // values are taken in increasing order
	      for (j=0; j<n; j++)
		  QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
	      
	      idxL--;idxL=FREE(idxL);
	      idxU--;idxU=FREE(idxU);

	      idxbuff =FREE(idxbuff);
	      idxposL =FREE(idxposL);
	      idxposU =FREE(idxposU);
	      invblock=FREE(invblock);
	      Ahead   =FREE(Ahead);
	      Alist   =FREE(Alist);
	      Afirst  =FREE(Afirst);
	      BLhead  =FREE(BLhead);
	      BLlist  =FREE(BLlist);
	      BLfirst =FREE(BLfirst);
	      BUThead =FREE(BUThead);
	      BUTlist =FREE(BUTlist);
	      BUTfirst=FREE(BUTfirst);
	      AvalD   =FREE(AvalD);
	      AvalL   =FREE(AvalL);
	      AvalU   =FREE(AvalU);
	      Arowind =FREE(Arowind);
	      Acolind =FREE(Acolind);
	      buffC   =FREE(buffC);
	      buff    =FREE(buff);

	      for (i=0; i<nblocks; i++) {
		  BL->colind[i] =FREE(BL->colind[i] );
		  BL->rowind[i] =FREE(BL->rowind[i] );
		  BUT->colind[i]=FREE(BUT->colind[i]);
		  BUT->rowind[i]=FREE(BUT->rowind[i]);
		  BiD->colind[i]=FREE(BiD->colind[i]);
		  BL->valE[i]   =FREE(BL->valE[i]   );
		  BUT->valE[i]  =FREE(BUT->valE[i]  );
		  BiD->valD[i]  =FREE(BiD->valD[i]  );
	      }
	      BL->nblockcol =FREE(BL->nblockcol);
	      BL->nblockrow =FREE(BL->nblockrow);
	      BL->colind    =FREE(BL->colind);
	      BL->rowind    =FREE(BL->rowind);
	      BL->valE      =FREE(BL->valE);

	      BUT->nblockcol=FREE(BUT->nblockcol);
	      BUT->nblockrow=FREE(BUT->nblockrow);
	      BUT->colind   =FREE(BUT->colind   );
	      BUT->rowind   =FREE(BUT->rowind   );
	      BUT->valE     =FREE(BUT->valE     );

	      BiD->nblockcol=FREE(BiD->nblockcol);
	      BiD->colind   =FREE(BiD->colind   );
	      BiD->valD     =FREE(BiD->valD     );
	      
	      BL->nc=BL->nr=BUT->nc=BUT->nr=BiD->nc=BiD->nr=0;
	      BL->nblocks=BUT->nblocks=BiD->nblocks=0;

	      if (flag_ilu1t) {
		 blocksize=FREE(blocksize);
		 blocksize=myblocksize;
		 nblocks  =mynblocks;
	      } // end if flag_ilu1t

	      if (perturbation)
		 buffDGLBLCK=FREE(buffDGLBLCK);
	      
	      return (ierr);
	   } // end if

	   
	   // compute local determinant
	   locdet.r=0.0;
	   locdet.i=0.0;
	   flag=0;
	   for (ii=0,jj=0; jj<bi; jj++,ii+=cntL) {
	       // toggle flag to compute the sign of the permutation
	       if (idxU[jj]!=jj+1)
		  flag=!flag;
	       // diagonal pivot
	       ajj=buff[ii+jj];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       // mexPrintf("ajj=%12.4le\n",ajj);
	       if (ajj>=0.0)
		  locdet.r+=log(ajj);
	       else { // log of negative numbers is complex (main branch)
		  locdet.r+=log(-ajj); 
		  locdet.i+=M_PI;
	       }
#else
	       // mexPrintf("ajj=%12.4le+i%12.4le\n",ajj.r,ajj.i);
	       // complex logarithm (main branch)
	       rval=FABS(ajj);
	       // 1. real part
	       locdet.r+=log(rval);
	       // 1. imaginary part
	       // compute argument of ajj
	       if (ajj.i>=0) 
		  rval=acos(ajj.r/rval);
	       else
		  rval=-acos(ajj.r/rval);
	       locdet.i+=rval;	       
#endif
	   } // end for jj
	   // sign change by permutation
	   if (flag)
	      locdet.i+=M_PI;
	   
	   if (invert_blocks) {
	      ii=1; // optimal block size for GETRI
	      jj=-1;
#if defined _SINGLE_REAL_
	      m=ILAENV(&ii, "SGETRI", transa, &bi, &jj, &jj, &jj,6,1);
#elif defined _DOUBLE_REAL_
	      m=ILAENV(&ii, "DGETRI", transa, &bi, &jj, &jj, &jj,6,1);
#elif defined _SINGLE_COMPLEX_
	      m=ILAENV(&ii, "CGETRI", transa, &bi, &jj, &jj, &jj,6,1);
#else
	      m=ILAENV(&ii, "ZGETRI", transa, &bi, &jj, &jj, &jj,6,1);
#endif
	      // use buffC as work array
	      l=m*bi;
	      if (l>nbuffC) {
		 nbuffC=ELBOW_BUFF*l+BUFF_EXT;
		 buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	      } // end if
	      GETRI(&bi, buff,&cntL, idxU, buffC,&l, &ierr);
	      if (ierr) {
		 if (ierr<0)
		    printf("LAPACK's GETRI: %ld-th argument had an illegal value\n", -ierr);
		 else
		    printf("LAPACK's GETRI routine encountered zero column in step %ld\n", ierr);

		 // sort original matrix back such that row indices and their
		 // values are taken in increasing order
		 for (j=0; j<n; j++)
		     QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
	
		 idxL--;idxL=FREE(idxL);
		 idxU--;idxU=FREE(idxU);

		 idxbuff =FREE(idxbuff);
		 idxposL =FREE(idxposL);
		 idxposU =FREE(idxposU);
		 invblock=FREE(invblock);
		 Ahead   =FREE(Ahead);
		 Alist   =FREE(Alist);
		 Afirst  =FREE(Afirst);
		 BLhead  =FREE(BLhead);
		 BLlist  =FREE(BLlist);
		 BLfirst =FREE(BLfirst);
		 BUThead =FREE(BUThead);
		 BUTlist =FREE(BUTlist);
		 BUTfirst=FREE(BUTfirst);
		 AvalD   =FREE(AvalD);
		 AvalL   =FREE(AvalL);
		 AvalU   =FREE(AvalU);
		 Arowind =FREE(Arowind);
		 Acolind =FREE(Acolind);
		 buffC   =FREE(buffC);
		 buff    =FREE(buff);

		 for (i=0; i<nblocks; i++) {
		     BL->colind[i] =FREE(BL->colind[i] );
		     BL->rowind[i] =FREE(BL->rowind[i] );
		     BUT->colind[i]=FREE(BUT->colind[i]);
		     BUT->rowind[i]=FREE(BUT->rowind[i]);
		     BiD->colind[i]=FREE(BiD->colind[i]);
		     BL->valE[i]   =FREE(BL->valE[i]   );
		     BUT->valE[i]  =FREE(BUT->valE[i]  );
		     BiD->valD[i]  =FREE(BiD->valD[i]  );
		 }
		 BL->nblockcol =FREE(BL->nblockcol);
		 BL->nblockrow =FREE(BL->nblockrow);
		 BL->colind    =FREE(BL->colind);
		 BL->rowind    =FREE(BL->rowind);
		 BL->valE      =FREE(BL->valE);

		 BUT->nblockcol=FREE(BUT->nblockcol);
		 BUT->nblockrow=FREE(BUT->nblockrow);
		 BUT->colind   =FREE(BUT->colind   );
		 BUT->rowind   =FREE(BUT->rowind   );
		 BUT->valE     =FREE(BUT->valE     );

		 BiD->nblockcol=FREE(BiD->nblockcol);
		 BiD->colind   =FREE(BiD->colind   );
		 BiD->valD     =FREE(BiD->valD     );
		 
		 BL->nc=BL->nr=BUT->nc=BUT->nr=BiD->nc=BiD->nr=0;
		 BL->nblocks=BUT->nblocks=BiD->nblocks=0;
	      
		 if (flag_ilu1t) {
		    blocksize=FREE(blocksize);
		    blocksize=myblocksize;
		    nblocks  =mynblocks;
		 } // end if flag_ilu1t

		 if (perturbation)
		    buffDGLBLCK=FREE(buffDGLBLCK);
	      
		 return (ierr);
	      } // end if 
	   } // end if invert_blocks
	   else {
	      // export permutation
	      memcpy(pivots,idxU,bi*sizeof(integer));
	   } // end if-else invert_blocks

	
	   // extract (inverse) diagonal block to BiD
	   BiD->valD[i_bldu]=(FLOAT *)malloc((size_t)bi*bi*sizeof(FLOAT));
	   pD=BiD->valD[i_bldu];
	   pb=buff;
	   for (ii=0; ii<bi; ii++,pD+=bi,pb+=cntL)
	       memcpy(pD,pb,(size_t)bi*sizeof(FLOAT));
#ifdef _PROFILING_
	   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // allocate memory for BL
	   ldc=cntL-bi; // number of rows before dropping
	   k=ldc*bi;    // total size
	   BL->valE[i_bldu]=(FLOAT *)malloc((size_t)k*sizeof(FLOAT));
	   if (invert_blocks) {
	      // multiply strict lower triangular block by the inverse diagonal block
	      // and copy it to BL
	      // use level-3-BLAS 
	      // init with 0
	      pL=BL->valE[i_bldu];
	      for (ii=0; ii<k; ii++,pL++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  *pL=0.0;
#else
		  pL->r=pL->i=0.0;
#endif
	      } // end for ii
	      // C = 1*A*B+0*C
	      transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      alpha=1.0;
	      beta =0.0;
#else
	      alpha.r=1.0; alpha.i=0.0;
	      beta.r =0.0; beta.i =0.0;
#endif
	      if (ldc && bi)
		 GEMM(transa,transb, &ldc,&bi,&bi, &alpha,
		      buff+bi,&cntL,
		      BiD->valD[i_bldu],&bi, &beta, BL->valE[i_bldu],&ldc,1,1);
	   } // end if invert_blocks
	   else { // LU case
	      // forward/backward solve using the LU factorization of the diagonal block
	      // copy subdiagonal part of buff to BL since solve works in-place with BL^T
	      if (k>nbuffC) {
		 nbuffC=ELBOW_BUFF*k+BUFF_EXT;
		 buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	      } // end if
	      // buffC = subdgl(buff)^T = subdgl([BiD(i,i);BL(:,i)])^T
	      // resp. buffC = subdgl(buff)^T = subdgl([BiD(i,i);BL(:,i)])^T
	      pb=buff+bi;
	      pL=buffC;
	      i2=1;
	      for (ii=0; ii<ldc; ii++,pb++,pL+=bi) 
		  COPY(&bi,pb,&cntL,pL,&i2);
	      // forward/backward solves (LU)^T * BL^T  <- BL^T
	      //                    <=>  (LU)^T * buffC <- buffC
	      GETRS("t", &bi, &ldc, BiD->valD[i_bldu],&bi, idxU, buffC,&bi, &ierr,1);
	      if (ierr<0) {
		 printf("LAPACK's GETRS: %ld-th argument had an illegal value\n", -ierr);
		 return (ierr);
	      }
	      // BL(:,i)=buffC^T
	      pb=buffC;
	      pL=BL->valE[i_bldu];
	      i2=1;
	      for (ii=0; ii<ldc; ii++,pb+=bi,pL++) 
		  COPY(&bi,pb,&i2,pL,&ldc);	      
	   } // end if-else invert_blocks

	   
	   // apply dropping to the strict lower triangular block of BL
	   pL=BL->valE[i_bldu];
	   pb=buff;   // remaining entries are buffered in "buff"
	   l=0;
	   for (ii=0; ii<ldc; ii++,pL++) {
	       jj=I_AMAX(&bi,pL,&ldc)-1;
	       k=jj*ldc;
	       // copy larger entries back
	       if (FABS(pL[k])>=droptol) {
		  COPY(&bi, pL,&ldc, pb,&ldc);
		  pb++;
		  idxU[l++]=idxL[ii+bi];
	       } // if
	   } // end for ii
	   BL->nblockrow[i_bldu]=l;
	   BL->rowind[i_bldu]=(integer *)malloc((size_t)l*sizeof(integer));
	   pL=BL->valE[i_bldu];
	   pb=buff;
	   pi=BL->rowind[i_bldu];
	   for (ii=0; ii<l; ii++,pL++,pb++) {
	       *pi++=idxU[ii];
	       COPY(&bi, pb,&ldc, pL,&l);
	   } // for ii
	   BL->valE[i_bldu]=(FLOAT *)realloc(BL->valE[i_bldu],(size_t)l*bi*sizeof(FLOAT));
#ifdef _PROFILING_
	   time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	} // end if-else bi=1

	
	// update determinant
	// mexPrintf("det=%12.4le+i%12.4le\n",locdet.r,locdet.i);
	determinant->r+=locdet.r;
	determinant->i+=locdet.i;

	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// reset idxposL
	for (m=0; m<cntL; m++) {
	    t=idxL[m];
	    idxposL[t]=0;
	} // end for m
	cntL=0;
#ifdef PRINT_CHECK
	for (m=0; m<n; m++) {
	  if (idxposL[m]) {
	    printf("8.idxposL[%ld]=%ld\n",m,idxposL[m]);
	    fflush(stdout);
	  } // end if
	  if (idxposU[m]) {
	    printf("8.idxposU[%ld]=%ld\n",m,idxposU[m]);
	    fflush(stdout);
	  } // end if
	} // end for m	
#endif
#ifdef _PROFILING_
	time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	// --- END compute block column i of BL using linked list of BUT ---
	// -----------------------------------------------------------------

	
	// now (almost) same procedure for BUT, already excluding the diagonal block

	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// -----------------------------------------------------------------
	// ----- compute block column i of BUT using linked list of BL -----
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 start
        r=BiD->colind[i_bldu][0]; s=r+bi;

	
	// first pass: set up nonzero row structure of U before dropping
	// only store indices of the super-diagonal block of A
	for (m=0; m<ncolind; m++) {
	    t=Acolind[m];
	    idxU[cntU]=t;
	    idxposU[t]=++cntU;
	} // end for m

	
	// beginning of the linked list for block column i
	// this information is obtained from the sparsity pattern of L
	it=BLhead[i_bldu];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning block column %ld of BUT\n",it); fflush(stdout);
#endif
	      // check block COLUMN "it" of BUT
	      rowind=BUT->rowind[it];
	      nrow=BUT->nblockrow[it];
	      // recall that we already advanced BUTfirst[it] such that its
	      // leading entry is now at least s
	      if (cntU) {
		 // idxU is non-empty, use idxL temporarily for the new indices
		 for (m=BUTfirst[it]; m<nrow; m++) {
		     t=rowind[m];
		     // skip entries above or in the diagonal block
		     if (t>=s) 
		        break;
		 } // end for m
		 for (; m<nrow; m++) {
		     t=rowind[m];
		     // only store new indices
		     if (!idxposU[t]) {
		        idxL[cntL++]=t;
		     } // new fill
		 } // end for m
		 // now idxU[0,...,cntU-1] and idxL[0,...,cntL-1] are
		 // distinct index arrays each with indices in increasing order
		 // we merge them in idxU
		 jj=cntU-1; ii=cntL-1;
		 i2=cntU+cntL;
		 tt=idxU[jj];
		 t=idxL[ii]; 
		 while (ii>=0) {
		       // new index greater than given index
		       if (t>tt) {
			  idxposU[t]=i2--;
			  idxU[i2]=t;
			  ii--;
			  t=idxL[ii]; 
		       }
		       else { // tt>t
			  while (tt>t) {
			        idxposU[tt]=i2--;
				idxU[i2]=tt;
				jj--;
				tt=idxU[jj];
			  } // end while
		       } // end if-else
		 } // end while
		 cntU+=cntL;
		 cntL=0;
	      }
	      else { // idxU has not been used yet, directly insert indices
		 // printf("cntL=%ld,cntU=%ld\n",cntL,cntU);fflush(stdout);
		 for (m=BUTfirst[it]; m<nrow; m++) {
		     t=rowind[m];
		     // skip entries above or in the diagonal block
		     if (t>=s)
		        break;
		 } // end for m
		 for (; m<nrow; m++) {
		     t=rowind[m];
		     // store new indices
		     idxU[cntU]=t;
		     idxposU[t]=++cntU;
		 } // end for m
	      } // end if-else cntU
		 
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		  if (idxposL[m]) {
		     printf("9.idxposL[%ld]=%ld\n",m,idxposL[m]);
		     fflush(stdout);
		  } // end if
	      } // end for m	
#endif
	      
	      // advance to next required block column of BUT as requested by BL
	      it=BLlist[it];
	} // end while
	// end first pass
#ifdef PRINT_INFO0
	printf("total nonzero row indices column %ld of BUT\n",i); fflush(stdout);
	for (m=0; m<cntU; m++) {
  	    printf("%8ld",idxU[m]);
	} // end for m
	printf("\n");fflush(stdout);
#endif
	
	// preparing second pass
	// now we know the size of the current block column before dropping
	// buff is going to store block column i of U temporarily
	if (cntU*bi>nbuff) {
	   // increase buffer
	   nbuff=ELBOW_BUFF*(cntU*bi)+BUFF_EXT;
	   buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	} // end if
	// init whole buff with 0 and scatter super-diagonal block of A into buff
 	pb=buff; pU=AvalU;
	for (m=0; m<bi; m++,pb+=cntU,pU+=ncolind) {
	    // init all rows with 0
	    for (jj=0; jj<cntU; jj++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	        pb[jj]=0.0;
#else
		pb[jj].r=pb[jj].i=0.0;
#endif
	    } // end for k
	    // copy columns of the super-diagonal block of A
	    for (jj=0; jj<ncolind; jj++) {
	        // column index t
	        t=Acolind[jj];
		// associated position shifted by 1
		tt=idxposU[t];
		pb[tt-1]=pU[jj];
	    } // end for jj
	} // end for m
#ifdef _PROFILING_
	time_LU_update_pass_1+=omp_get_wtime()-timeBegin;
#endif
	
	// second pass: compute block column i of BUT using level-3-BLAS
	// beginning of the linked list of BL, the nonzero pattern of BL
	// reflects the requested columns of BUT
	it=BLhead[i_bldu];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning again block column %ld of BUT\n",it); fflush(stdout);
#endif
	      // level-3-BLAS update C+=AB^T

#ifdef _PROFILING_
	      timeBegin = omp_get_wtime();
#endif
	      // step 1: rows of block column "it" of BL refer to rows of B
	      //         These can be left in-place for the level-3-BLAS update
	      // read the associated rows from BL for B
	      jj=BLfirst[it];        // start of row r and higher
	      ldb=BL->nblockrow[it]; // leading dimension (total number of rows)
	      flag=-1;
	      rowind=BL->rowind[it]; // shortcut
	      while (jj<ldb && flag) {
		    // row indices between r and s-1
		    t=rowind[jj];
		    // row indices are inside the current block column i
		    if (t<s) {
		       idxL[cntL]=t;      // bookmark index and ...
		       idxposL[t]=++cntL; // its position
		       jj++;
		    }
		    else // block column i exceeded
		       flag=0;
	      } // end while
	      // now idxL, idxposL refer to the columns of block column i of BUT
	      // to be updated with "it"
#ifdef PRINT_INFO0
	      printf("column indices required to update with column %ld of BUT\n",it); fflush(stdout);
	      for (m=0; m<cntL; m++) {
		  printf("%8ld",idxL[m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef PRINT_INFO0
	      printf("row indices required to update with column %ld of BUT\n",it); fflush(stdout);
	      for (m=BUTfirst[it]; m<BUT->nblockrow[it]; m++) {
		  printf("%8ld",BUT->rowind[it][m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	      time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif

	      l=BUT->nblockcol[it];
	      ii=BUTfirst[it];        // start of row s and higher in BUT
	      lda=BUT->nblockrow[it]; // leading dimension of block column "it" of BUT
	      ldc=lda-ii;             // number of rows to be updated
	      if (ldc && cntL && l) {
#ifdef _PROFILING_
		 timeBegin=omp_get_wtime();
#endif
		 // single-column-update, do this by hand and skip level-3-BLAS
		 // compress index arrays into one array
		 pi=BUT->rowind[it]+ii; // row indices to be considered
		 for (k=0; k<ldc; k++) {
		     // which row indices are required?
		     tt=*pi++;
		     idxbuff[k]=idxposU[tt];
		 } // end for k
#ifdef _PROFILING_
		 time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif
		 if (cntL==1 && l==1) {
#ifdef _PROFILING_
		    timeBegin=omp_get_wtime();
#endif
		    // extract column index t
		    t=idxL[0];
		    // scatter to column t of buff
		    pb=buff+cntU*(t-r)-1;   // associated column of buff,shifted
		    pU=BUT->valE[it]  +ii; // numerical values of column BUT{it}
		    pL=BL->valE[it]   +BLfirst[it]; // numerical value BL{it}
		    
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
		    val=-(*pL);
#else
		    val.r=-pL->r;
		    val.i=-pL->i;
#endif

#ifdef _USE_MKL_
		    // use sparse AXPYI pb[idxbuff]= val * pU + pb[idxbuff]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    CBLAS_AXPYI(ldc, val,  pU,idxbuff, pb); 
#else
		    CBLAS_AXPYI(ldc, &val, pU,idxbuff, pb); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
		    
#else // hand-coded loop
		    for (k=0; k<ldc; k++) {
			// position of tt in buff, shifted by 1
			i2=idxbuff[k];
			// downdate buff - BUT{i}*BL{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			pb[i2]+=(*pU++)*val;
#else
			pb[i2].r+=pU->r*val.r-pU->i*val.i;
			pb[i2].i+=pU->r*val.i+pU->i*val.r;
			pU++;
#endif
		    } // end for k
#endif //-else _USE_MKL_
#ifdef _PROFILING_
		    time_scalar_update+=omp_get_wtime()-timeBegin;
#endif
		    
		 }
		 else if ((double)ldc*cntL<=SMALL_BLOCK &&
			  (double)ldc*l   <=SMALL_BLOCK &&
			  (double)l*cntL  <=SMALL_BLOCK) {
 		    // use hand-coded loops and update in-place

		    // rank-1-update
		    pL=BL->valE[it] +BLfirst[it]; // numerical value BL{it}
		    if (l==1) {
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
 		       for (m=0; m<cntL; m++,pL++) {
			   // numerical values of column BUT{it}
			   pU=BUT->valE[it]  +ii; 
			   // extract column index t
			   t=idxL[m];
			   // directly scatter to column t of buff
			   pb=buff+cntU*(t-r)-1;   // associated column of buff, shifted

#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   val=-(*pL);
#else
			   val.r=-pL->r;
			   val.i=-pL->i;
#endif
			   
#ifdef _USE_MKL_
			   // use sparse AXPYI pb[idxbuff]= val * pU + pb[idxbuff]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   CBLAS_AXPYI(ldc, val,  pU,idxbuff, pb); 
#else
			   CBLAS_AXPYI(ldc, &val, pU,idxbuff, pb); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
		    
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       // downdate buff - BUT{i}*BL{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			       pb[i2]+=(*pU++)*val;
#else
			       pb[i2].r+=pU->r*val.r-pU->i*val.i;
			       pb[i2].i+=pU->r*val.i+pU->i*val.r;
			       pU++;
#endif
			   } // end for k
#endif //-else _USE_MKL_

		       } // end for m
#ifdef _PROFILING_
		       time_rank_1_update+=omp_get_wtime()-timeBegin;
#endif
		    }
		    else { // higher rank using DOT product
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
		       for (m=0; m<cntL; m++,pL++) {
			   // numerical values of column BUT{it}
			   pU=BUT->valE[it]  +ii; 
			   // extract column index t
			   t=idxL[m];
			   // directly scatter to column t of buff
			   pb=buff+cntU*(t-r)-1; // associated column of buff, shifted
			   
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       // downdate buff - BUT{i}*BL{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			       pb[i2]-=DISTDOT(&l, pU,&lda, pL,&ldb);
#else
#ifdef _USE_MKL_
			       DISTDOTU(&val, &l, pU,&lda, pL,&ldb);
#else
			       val=DISTDOTU(&l, pU,&lda, pL,&ldb);
#endif
			       pb[i2].r-=val.r;
			       pb[i2].i-=val.i;
#endif
			       pU++;
			   } // end for k
		       } // end for m
#ifdef _PROFILING_
		       time_small_block_update+=omp_get_wtime()-timeBegin;
#endif
		    } // end if-else l=1
		 }
		 else {
#ifdef _PROFILING_
		    timeBegin=omp_get_wtime();
#endif
		    // check whether we are able to bypass the gather/scatter part
		    // this is the case, whenever the selected column indices are
		    // a contiguous sequence and the associated row indices are
		    // contiguous sequence as well
		    if (idxbuff[ldc-1]-idxbuff[0]+1==ldc &&
			idxL[cntL-1]  -idxL[0]   +1==cntL) {
		       skip_gthr_sctr=-1;
		       // printf("gather/scatter could be skipped\n"); fflush(stdout);
		    }
		    else
		       skip_gthr_sctr=0;

		    // step 2: temporarily gather the data from buff to (buff)C to
		    //         accomodate the level-3-BLAS update
		    if (!skip_gthr_sctr) {
		       // extract the rows of buff associated with the rows of BUT
		       // recall again that we already excluded the diagonal block by
		       // advancing BUTfirst[it]
		       if (cntL*ldc>nbuffC) {
			  nbuffC=ELBOW_BUFF*(cntL*ldc)+BUFF_EXT;
			  buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
		       } // end if
	      
		       // gather the associated rows and columns from buff to (buff)C
		       pC=buffC;
		       for (m=0; m<cntL; m++) {
			   // extract column index t
			   t=idxL[m];
			   // gather column t from buff into buffC
			   pb=buff+cntU*(t-r)-1; // associated column of buff, shifted
			
#ifdef _USE_MKL_
			   CBLAS_GTHR(ldc, pb, pC,idxbuff);
			   pC+=ldc;
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       *pC++=pb[i2];
			   } // end for k
#endif //-else _USE_MKL_
			
		       } // end for m
		    } // end if NOT skip_gthr_sctr
#ifdef _PROFILING_
		    time_level_3_blas_gather+=omp_get_wtime()-timeBegin;
#endif

		    // now perform level-3-BLAS update
		    // C = -1*A*B^T+1*C
		    transa="n"; transb="t";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    alpha=-1.0;
		    beta = 1.0;
#else
		    alpha.r=-1.0; alpha.i=0.0;
		    beta.r = 1.0; beta.i =0.0;
#endif
		    // level-3 BLAS performed in-place
		    if (skip_gthr_sctr)
		       GEMM(transa,transb, &ldc,&cntL,&l, &alpha,
			    BUT->valE[it]+BUTfirst[it],&lda,
			    BL->valE[it]+BLfirst[it],&ldb, &beta,
			    buff + idxbuff[0]-1 + cntU*(idxL[0]-r),&cntU, 1,1);
		    else // level-3 BLAS uses a cached copy
		       GEMM(transa,transb, &ldc,&cntL,&l, &alpha,
			    BUT->valE[it]+BUTfirst[it],&lda,
			    BL->valE[it]+BLfirst[it],&ldb, &beta,
			    buffC,&ldc, 1,1);
	      

#ifdef _PROFILING_
		    timeBegin2=omp_get_wtime();
#endif
		    if (!skip_gthr_sctr) {
		       // scatter the result back to the associated rows
		       // and columns of buff
		       pC=buffC;
		       for (m=0; m<cntL; m++) {
		           // extract column index t
			   t=idxL[m];
			   // scatter to column t of buff
			   pb=buff+cntU*(t-r)-1; // associated column of buff, shifted

#ifdef _USE_MKL_
			   CBLAS_SCTR(ldc, pC, idxbuff, pb);
			   pC+=ldc;
#else // hand-coded loop
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       pb[i2]=*pC++;
			   } // end for k
#endif //-else _USE_MKL_

		       } // end for m
		    } // end if NOT skip_gthr_sctr
#ifdef _PROFILING_
		    time_level_3_blas_update+=omp_get_wtime()-timeBegin;
		    time_level_3_blas_scatter+=omp_get_wtime()-timeBegin2;
#endif
		 } // end if cntL=1 and l=1
	      } // end if ldc>0 and cntL>0 and l>0

#ifdef _PROFILING_
	      timeBegin=omp_get_wtime();
#endif
	      // update BLfirst, this automatically advances the beginning
	      // to the first position behind the diagonal block
	      BLfirst[it]=jj;      

	      // reset idxposL
	      for (m=0; m<cntL; m++) {
		  t=idxL[m];
		  idxposL[t]=0;
	      } // end for m
	      cntL=0;
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		  if (idxposL[m]) {
		     printf("10.idxposL[%ld]=%ld\n",m,idxposL[m]);
		     fflush(stdout);
		  } // end if
	      } // end for m	
#endif
	      
	      // advance to the next block column
	      it=BLlist[it];
#ifdef _PROFILING_
	      time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	} // end while
	// end second pass
  
	
	
	// scalar column, skip LAPACK
	if (bi==1) {
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // apply dropping to the strict lower triangular block of BUT
	   // cntU = number of rows before dropping
	   // rval = |diagonal_entry|*droptol
	   pb=buff; // start strict lower triangular part
	   l=0;
	   for (ii=0; ii<cntU; ii++,pb++)
	       if (FABS(*pb)>=rval)
		  l++;
	   // now we know the fill
	   BUT->nblockrow[i_bldu]=l;
	   BUT->valE[i_bldu]=(FLOAT *)malloc((size_t)l*sizeof(FLOAT));
	   BUT->rowind[i_bldu]=(integer *)malloc((size_t)l*sizeof(integer));
	   	   
	   // copy remaining strict lower triangular part to BUT
	   pb=buff;
	   pU=BUT->valE[i_bldu];
	   pi=BUT->rowind[i_bldu];
	   for (ii=0; ii<cntU; ii++,pb++)
	       if (FABS(*pb)>=rval) {
		  *pi++=idxU[ii];
		  *pU++=*pb;
	       }
#ifdef _PROFILING_
	       time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	}
	else { // bi>1, block case
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // allocate memory for BUT
	   k=cntU*bi;    // total size
	   BUT->valE[i_bldu]=(FLOAT *)malloc((size_t)k*sizeof(FLOAT));
	   if (invert_blocks) {
	      // multiply (transposed) strict upper triangular block by the inverse
	      // diagonal block and copy it to BUT
	      // use level-3-BLAS 
	      // init with 0
	      pU=BUT->valE[i_bldu];
	      for (ii=0; ii<k; ii++,pU++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	          *pU=0.0;
#else
		  pU->r=pU->i=0.0;
#endif
	      } // end for ii
	      // C = 1*A*B^T+0*C
	      transa="n"; transb="t";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      alpha=1.0;
	      beta =0.0;
#else
	      alpha.r=1.0; alpha.i=0.0;
	      beta.r =0.0; beta.i =0.0;
#endif
	      if (cntU && bi)
		 GEMM(transa,transb, &cntU,&bi,&bi, &alpha,
		      buff,&cntU,
		      BiD->valD[i_bldu],&bi, &beta, BUT->valE[i_bldu],&cntU,1,1);
	   } // end if invert_blocks
	   else { // LU case
	      // forward/backward solve using the LU factorization of the diagonal block
	      // copy subdiagonal part of buff to BUT since solve works in-place with BUT^T
	      if (k>nbuffC) {
		 nbuffC=ELBOW_BUFF*k+BUFF_EXT;
		 buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	      } // end if
	      // buffC = BUT(:,i)^T
	      // resp. buffC = buff^T 
	      pb=buff;
	      pU=buffC;
	      i2=1;
	      for (ii=0; ii<cntU; ii++,pb++,pU+=bi)
		  // pb -> pU
		  COPY(&bi,pb,&cntU,pU,&i2);
	      // forward/backward solves LU * BUT^T  <- BUT^T
	      //                    <=>  LU * buffC <- buffC
	      GETRS("n", &bi, &cntU, BiD->valD[i_bldu],&bi, pivots, buffC,&bi, &ierr,1);
	      if (ierr<0) {
		 printf("LAPACK's GETRS: %ld-th argument had an illegal value\n", -ierr);
		 return (ierr);
	      }
	      // BUT(:,i)=buffC^T
	      pb=buffC;
	      pU=BUT->valE[i_bldu];
	      i2=1;
	      for (ii=0; ii<cntU; ii++,pb+=bi,pU++) 
		  // pb -> pU
		  COPY(&bi,pb,&i2,pU,&cntU);	      
	   } // end if-else invert_blocks

	   
	   // apply dropping to the (transposed) strict upper triangular block of BUT
	   pU=BUT->valE[i_bldu];
	   l=0;
	   for (ii=0; ii<cntU; ii++,pU++) {
	       jj=I_AMAX(&bi,pU,&cntU)-1;
	       k=jj*cntU;
	       // copy indices of larger entries
	       if (FABS(pU[k])>=droptol) {
		  idxL[l++]=idxU[ii];
	       } // if
	   } // end for ii
	   BUT->nblockrow[i_bldu]=l;
	   BUT->rowind[i_bldu]=(integer *)malloc((size_t)l*sizeof(integer));
	   pU=BUT->valE[i_bldu];
	   pi=BUT->rowind[i_bldu];
	   for (ii=0; ii<l; ii++,pU++) {
	       t=idxL[ii]; // ==idxU[for some index m]
	       m=idxposU[t]-1; // position shifted by 1
	       *pi++=t;
	       // now copy original U entries without inverse to BUT!
	       COPY(&bi, buff+m,&cntU, pU,&l);
	   } // for ii
	   BUT->valE[i_bldu]=(FLOAT *)realloc(BUT->valE[i_bldu],(size_t)l*bi*sizeof(FLOAT));
#ifdef _PROFILING_
	   time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	} // end if-else bi=1
	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// reset idxposU
	for (m=0; m<cntU; m++) {
	    t=idxU[m];
	    idxposU[t]=0;
	} // end for m
	cntU=0;
#ifdef PRINT_CHECK
	for (m=0; m<n; m++) {
	  if (idxposU[m]) {
	    printf("11.idxposU[%ld]=%ld\n",m,idxposU[m]);
	    fflush(stdout);
	  } // end if
	} // end for m	
#endif
#ifdef _PROFILING_
	time_LU_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	// --- END compute block column i of BUT using linked list of BL ---
	// -----------------------------------------------------------------
#ifdef PRINT_INFO0
	printf("block %3ld of BL, BiD, BUT\n",i);
	if (invert_blocks)
	   printf("inverse diagonal block BiD[%3ld]\n        ",i);
	else
	   printf("factorized diagonal block BiD[%3ld]\n        ",i);
	for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	    printf("%8ld",BiD->colind[i_bldu][m]);
	printf("\n");
	for (k=0; k<BiD->nblockcol[i_bldu]; k++) {
	    printf("%8ld",BiD->colind[i_bldu][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",BiD->valD[i_bldu][m*BiD->nblockcol[i_bldu]+k]);
	    printf("\n");
#else	    
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",BiD->valD[i_bldu][m*BiD->nblockcol[i_bldu]+k].r);
	    printf("\n        ");
	    for (m=0; m<BiD->nblockcol[i_bldu]; m++) 
	        printf("%8.1le",BiD->valD[i_bldu][m*BiD->nblockcol[i_bldu]+k].i);
	    printf("\n");
#endif
	}
	printf("\n");fflush(stdout);
	if (!invert_blocks) {
	   printf("pivots: ");
	   for (m=0; m<BiD->nblockcol[i_bldu]; m++)
	       printf("%8ld",pivots[m]);
	   printf("\n");fflush(stdout);
	}

	printf("sub-diagonal block BL\n        ");
	if (BL->nblockrow[i_bldu]>0) {
	   for (m=0; m<BL->nblockcol[i_bldu]; m++) 
	       printf("%8ld",BL->colind[i_bldu][m]);
	   printf("\n");
	   for (k=0; k<BL->nblockrow[i_bldu]; k++) {
	       printf("%8ld",BL->rowind[i_bldu][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       for (m=0; m<BL->nblockcol[i_bldu]; m++) 
	           printf("%8.1le",BL->valE[i_bldu][m*BL->nblockrow[i_bldu]+k]);
#else
	       for (m=0; m<BL->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",BL->valE[i_bldu][m*BL->nblockrow[i_bldu]+k].r);
	       printf("\n        ");
	       for (m=0; m<BL->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",BL->valE[i_bldu][m*BL->nblockrow[i_bldu]+k].i);
#endif
	       printf("\n");
	   } // end for k
	}
	printf("\n");

	printf("transposed super-diagonal block BUT\n        ");
	if (BUT->nblockrow[i_bldu]>0) {
	   for (m=0; m<BUT->nblockcol[i_bldu]; m++) 
	       printf("%8ld",BUT->colind[i_bldu][m]);
	   printf("\n");
	   for (k=0; k<BUT->nblockrow[i_bldu]; k++) {
	       printf("%8ld",BUT->rowind[i_bldu][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       for (m=0; m<BUT->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",BUT->valE[i_bldu][m*BUT->nblockrow[i_bldu]+k]);
#else
	       for (m=0; m<BUT->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",BUT->valE[i_bldu][m*BUT->nblockrow[i_bldu]+k].r);
	       printf("\n        ");
	       for (m=0; m<BUT->nblockcol[i_bldu]; m++) 
		   printf("%8.1le",BUT->valE[i_bldu][m*BUT->nblockrow[i_bldu]+k].i);
#endif
	       printf("\n");
	   } // end for k
	}
	printf("\n");
#endif


#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// check fill-in of current and previous block in order to find out
	// whether a progressive aggregation is recommended
	if (i_bldu && PROGRESSIVE_AGGREGATION && pa) {
	   // first compute old fill
	   k=BiD->nblockcol[i_bldu-1];
	   l=BL->nblockrow[i_bldu-1];
	   m=BUT->nblockrow[i_bldu-1];
	   ii=k*(k+l+m);
	   r=BiD->nblockcol[i_bldu];
	   s=BL->nblockrow[i_bldu];
	   t=BUT->nblockrow[i_bldu];
	   ii+=r*(r+s+t);

	   // now determine new fill
	   // checkmark fill of the current sub-diagonal block
	   pi=BL->rowind[i_bldu];
	   for (jj=0; jj<s; jj++) {
	       idxL[cntL]=*pi;
	       idxposL[*pi++]=++cntL;
	   } // end for jj
	   // compute additional fill from the previous sub-diagonal block
	   jj=BLfirst[i_bldu-1];
	   pi=BL->rowind[i_bldu-1]+jj;
	   for (; jj<l; jj++) {
	       // only insert additional fill
	       if (!idxposL[*pi]) {
		  idxL[cntL]=*pi;
		  idxposL[*pi]=++cntL;
	       }
	       pi++;
	   } // end for jj
	   // checkmark fill of the current transposed super-diagonal block
	   pi=BUT->rowind[i_bldu];
	   for (jj=0; jj<t; jj++) {
	       idxU[cntU]=*pi;
	       idxposU[*pi++]=++cntU;
	   } // end for jj
	   // compute additional fill from the previous transposed super-diagonal block
	   jj=BUTfirst[i_bldu-1];
	   pi=BUT->rowind[i_bldu-1]+jj;
	   for (; jj<m; jj++) {
	       // only insert additional fill
	       if (!idxposU[*pi]) {
		  idxU[cntU]=*pi;
		  idxposU[*pi]=++cntU;
	       }
	       pi++;
	   } // end for jj
	   
	   // now we are able to compute the new fill
	   jj=(k+r)*(k+r+cntL+cntU);



	   // we may merge these two block columns without too much additional fill
	   // we allow an additional percentage fill or
	   // a small additional number of rows
	   // we also claim that the new diagonal block size must 
	   // 1) be beyond a lower bound
	   // 2) be below an upper bound
	   // 3) grow substantially
	   if ((     jj<=ELBOW_BLOCK*(ii-cnt_fill)
	          || jj<=(ii-cnt_fill)+BLOCK_EXT*MAX(k,r))
	       && k+r>=MIN_BLOCK_SIZE_JANUS && k+r<=MAX_BLOCK_SIZE_JANUS
	       && k+r>=GROW_BLOCK_SIZE_JANUS*MAX(k,r)) {
	      // printf("block column %2ld: old fill %8ld, new fill %8ld\n",i_bldu,ii-cnt_fill,jj);

	      // incorporate additional fill-in
	      cnt_fill+=jj-ii;
	      // leading dimension of the enlarged BiD->valD[i_bldu]
              ldb=k+r;

	      // enlarge blocks
	      // 1. BiD
	      // BiD will now have some diagonal block of type
	      // [D11_new D12_new]
	      // [D21_new D22    ]
	      // with blocks to be computed
	      BiD->colind[i_bldu]=(integer *)realloc(BiD->colind[i_bldu],
						     (size_t)ldb*sizeof(integer));
	      BiD->valD[i_bldu]  =(FLOAT *)  realloc(BiD->valD[i_bldu],
						     (size_t)ldb*ldb*sizeof(FLOAT));

	      // 2. BL
	      // BL->colind[i_bldu-1] and BL->valE[i_bldu-1] are interpreted
	      // as some kind of block matrix
	      // [L21]
	      // [L31]
	      // where L21 is associated with block i_bldu (i.e. D22)
	      // and L31 refers to the bottom part
	      // BL->colind[i_bldu] and BL->valE[i_bldu] are now extended to cover
	      // a bigger sub-diagonal block of type
	      // [L31_new  L32]
	      // where L31_new is to be computed from the old block structures
	      // and L32 might be interlaced with a few zero rows
	      BL->colind[i_bldu]=(integer *)realloc(BL->colind[i_bldu],
						    (size_t)ldb*sizeof(integer));
	      BL->rowind[i_bldu]=(integer *)realloc(BL->rowind[i_bldu],
						    (size_t)cntL*sizeof(integer));
	      BL->valE[i_bldu]  =(FLOAT *)  realloc(BL->valE[i_bldu],
						    (size_t)cntL*ldb*sizeof(FLOAT));

	      // 3. BUT
	      // BUT->colind[i_bldu-1] and BUT->valE[i_bldu-1] are interpreted
	      // as some kind of block matrix
	      // [BUT21]  = [U12 U13]^T 
	      // [BUT31]
	      // where U11 is associated with block i_bldu (i.e. D22)
	      // and U13 refers to the bottom part
	      // BUT->colind[i_bldu] and BUT->valE[i_bldu] are now extended to cover
	      // a bigger sub-diagonal block of type
	      // [BUT31_new  BUT32] = [U13_new]^T
	      //                      [U23    ]
	      // where U13_new is to be computed from the old block structures
	      // and BUT32 might be interlaced with a few zero rows
	      BUT->colind[i_bldu]=(integer *)realloc(BUT->colind[i_bldu],
						     (size_t)ldb*sizeof(integer));
	      BUT->rowind[i_bldu]=(integer *)realloc(BUT->rowind[i_bldu],
						     (size_t)cntU*sizeof(integer));
	      BUT->valE[i_bldu]  =(FLOAT *)  realloc(BUT->valE[i_bldu],
						     (size_t)cntU*ldb*sizeof(FLOAT));

	      
	      
	      // we have to keep in mind that at the moment, BL==L
	      // is block unit lower triangular but BUT=U^T is not, the
	      // multiplication with D11, D22 is still missing (for
	      // performance reasons to faster downdate BL,BUT), i.e. we have
	      // a local factorization of type
	      // in the case of invert_blocks:
	      //   [ I   0 ] [D11^{-1}   U12    U13]   
	      //   [L21  I ] [  0      D22^{-1} U23] 
	      //   [L31 L32]                          
	      //
	      //   [ I   0 ] [D11^{-1}   0     ] [I D11*BUT12^T D11*BUT13^T]
	      // = [L21  I ] [  0      D22^{-1}] [0      I      D22*BUT23^T]
	      //   [L31 L32]                        
	      // respectively in the other case we have
	      //   [ I   0 ] [P11  0 ] [L11  0 ] [U11  0 ] [I (P11L11U11)^{-1}*BUT12^T (P11L11U11)^{-1}*BUT13^T]
	      // = [L21  I ] [ 0  P22] [ 0  L22] [ 0  U22] [0               I          (P22L22U22)^{-1}*BUT23^T]
	      //   [L31 L32]                        
	      

	      
	      // ----------------------------------------------------
	      // compute aggregated lower triangular block
	      // ----------------------------------------------------

	      // [   I     0 ]   [     I       0 ]   [ I   0 ] [  I  0]
	      // [   0     I ] = [     0       I ] = [L21  I ] [-L21 I]
	      // [L31_new L32]	 [L31-L32*L21 L32]   [L31 L32] 
	      
	      // 1. shift and copy

	      // shift L32 and supplement it with additional zero rows
	      // position of the last entry of L32 in the old block
	      pL=BL->valE[i_bldu]+s*r-1;
	      // future position of the last entry of the extended L32 inside
	      // the enlarged block
	      pb=BL->valE[i_bldu]+cntL*ldb-1;
	      // shift column-by-column
	      for (i2=0; i2<r; i2++) {
	          // position of the last fill row inside idxL
		  ii=cntL-1; 
		  // position of the last existing row inside idxL
		  jj=s-1;
		  while (ii>=s && jj>=0) {
	                // currently greatest fill index
			mm=idxL[ii];
			// currently greatest existing index
			tt=idxL[jj];
			// insert 0
			if (mm>tt) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                           *pb=0.0;
#else
                           pb->r=pb->i=0.0;
#endif		      
  		           ii--;
			}
			// we must have mm!=tt by construction! 
			else { // mm<tt, shift existing entry
		           *pb=*pL--;
                           jj--;
		        } // end if-else
			pb--;
	          } // end while
		  while (ii>=s) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                        *pb=0.0;
#else
                        pb->r=pb->i=0.0;
#endif		      
  		        ii--;
			pb--;
	          } // end while
		  while (jj>=0) {
		        *pb--=*pL--;
			jj--;
	          } // end while
	      } // end for i2

	      // copy L31 from the previous block and supplement it with zero rows
	      // position of the last entry of L31
	      pL=BL->valE[i_bldu-1]+l*k-1;
	      // future position of the last entry of the extended L31 inside
	      // the enlarged block
	      pb=BL->valE[i_bldu]+cntL*k-1;
	      // index array column i_bldu-1
	      pi=BL->rowind[i_bldu-1];
	      // index array column i_bldu
	      pi2=BL->rowind[i_bldu];
	      // number of rows of L21
	      it=BLfirst[i_bldu-1];
	      for (i2=0; i2<k; i2++,pL-=it) {
	          // position of the last row inside BL->nblockrow[i_bldu-1]
		  ii=l-1; 
	          // position of the last row inside BL->nblockrow[i_bldu]
		  jj=s-1;
		  while (ii>=it && jj>=0) {
	                // currently largest fill index
			mm=pi[ii];
			// currently largest existing index
			tt=pi2[jj];
			// insert 0
			if (tt>mm) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                           *pb=0.0;
#else
                           pb->r=pb->i=0.0;
#endif		      
  		           jj--;
			}
			else if (tt<mm) { // copy existing entry
		           *pb=*pL--;
                           ii--;
		        } 
			else { // common index, copy and decrement both counters
		           *pb=*pL--;
                           ii--;
                           jj--;
		        } // end if-else
			pb--;
	          } // end while
		  while (ii>=it) {
		        *pb--=*pL--;
                        ii--;
	          } // end while
		  while (jj>=0) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                        *pb=0.0;
#else
                        pb->r=pb->i=0.0;
#endif		      
  		        jj--;
			pb--;
	          } // end while
	      } // end for i2
	      
	      // 2. update L31_new = L31-L32*L21
	      // truncate columns of L32 to those associated with the rows of L21
	      // beginning of L32 within the enlarged block
	      pL=BL->valE[i_bldu]+cntL*k;
	      // array of nonzero indices of L21
	      pi=BL->rowind[i_bldu-1];
	      // first index of the diagonal block, we need this information in
	      // order to correctly adjust the pointer to the columns of L32
	      tt=BL->colind[i_bldu][0];
	      // number of columns to be cached
	      mm=BLfirst[i_bldu-1];
	      // use buff to hold the copy
	      if (cntL*m>nbuff) {
		 // increase buffer
		 nbuff=ELBOW_BUFF*(cntL*m)+BUFF_EXT;
		 buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	      } // end if
	      pb=buff;
	      i2=1;
	      for (ii=0; ii<mm; ii++,pi++,pb+=cntL)
	          // copy logical column *pi of L32 to buff
	          COPY(&cntL, pL+(*pi-tt)*cntL,&i2, pb,&i2);	          
	      // now we are ready to apply GEMM:
	      // L31_new =   -L32 *L21+  L31
	      //         = -1*buff*L21+1*L31_new
	      transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
              alpha=-1.0;   beta=1.0;
#else
              alpha.r=-1.0; beta.r=1.0;
	      alpha.i= 0.0; beta.i=0.0;
#endif
              if (cntL && k && mm)
                 GEMM(transa,transb, &cntL,&k,&mm, &alpha,
                      buff,&cntL,
		      BL->valE[i_bldu-1],&l, &beta,
		      BL->valE[i_bldu],&cntL, 1,1);


	      // 3. finalization
	      // to finalize the aggregated sub-diagonal block, we extend
	      // the column index set as well as the row index set
	      // shift existing indices
	      // old index position last column index
	      pi=BL->colind[i_bldu]+r-1;
	      // new index position last column index
	      pi2=BL->colind[i_bldu]+ldb-1;
	      for (ii=0; ii<r; ii++)
	          *pi2--=*pi--;
	      // insert column indices from the previous block
	      memcpy(BL->colind[i_bldu],BL->colind[i_bldu-1],k*sizeof(integer));
	      BL->nblockcol[i_bldu]=ldb;
	      // insert new row indices in increasing order
	      // position of the last fill row inside idxL
	      ii=cntL-1; 
	      // position of the last existing row inside idxL
	      jj=s-1;
	      // future position of the last index of the extended index set inside
	      // the enlarged block
	      pi=BL->rowind[i_bldu]+cntL-1;
	      while (ii>=s && jj>=0) {
	            // currently largest fill index
	            mm=idxL[ii];
		    // currently largest existing index
		    tt=idxL[jj];
		    // insert mm
		    if (mm>tt) {
		       *pi=mm;
		       ii--;
		    }
		    else { // shift existing index (here copied back from idxL)
		       *pi=tt;
                       jj--;
		    } // end if-else
		    pi--;
	      } // end while
	      while (ii>=s) {
	            // currently largest fill index
	            mm=idxL[ii];
		    *pi--=mm;
		    ii--;
	      } // end while
	      while (jj>=0) {
		    // currently largest existing index
		    tt=idxL[jj];
		    *pi--=tt;
		    jj--;
	      } // end while
	      BL->nblockrow[i_bldu]=cntL;

#ifdef PRINT_INFO0
	      printf("extended BL{%ld}\n",i_bldu);
	      if (BL->nblockrow[i_bldu]) {
		 printf("        ");
		 for (jj=0; jj<BL->nblockcol[i_bldu]; jj++) {
		     printf("%8ld",BL->colind[i_bldu][jj]);
		 }
		 printf("\n");
		 for (ii=0; ii<BL->nblockrow[i_bldu]; ii++) {
		     printf("%8ld",BL->rowind[i_bldu][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     for (jj=0; jj<BL->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BL->valE[i_bldu][BL->nblockrow[i_bldu]*jj+ii]);
		     }
#else
		     for (jj=0; jj<BL->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BL->valE[i_bldu][BL->nblockrow[i_bldu]*jj+ii].r);
		     }
		     printf("\n        ");		  
		     for (jj=0; jj<BL->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BL->valE[i_bldu][BL->nblockrow[i_bldu]*jj+ii].r);
		     }
#endif
		     printf("\n");		  
		 }
		 printf("\n");
	      }
#endif	      
	      // ----------------------------------------------------
	      // END compute aggregated lower triangular block
	      // ----------------------------------------------------


	      
	      // --------------------------------------------------------
	      // compute aggregated transposed upper triangular block
	      // --------------------------------------------------------

	      // [  I        0    ]   [  I            0        ]
	      // [  0        I    ] = [  0            I        ]
	      // [BUT31  BUT32_new]   [BUT31  BUT32+BUT31*L21^T]
	      // This is because we have
	      //   [ [D11_new D12_new]{-T}]
	      //   [ [D21_new D22_new]    ]
	      //   [BUT31        BUT32_new]
	      //
	      //   [D11^{-T}    0    ] [I L21^T]
	      // = [ BUT21   D22^{-T}] [0   I  ]
	      //   [ BUT31    BUT32  ] 	       
	      
	      // 1. shift and copy

	      // shift BUT32 and supplement it with additional zero rows
	      // position of the last entry of BUT32 in the old block
	      pU=BUT->valE[i_bldu]+t*r-1;
	      // future position of the last entry of the extended BUT32 inside
	      // the enlarged block
	      pb=BUT->valE[i_bldu]+cntU*ldb-1;
	      // shift column-by-column
	      for (i2=0; i2<r; i2++) {
	          // position of the last fill row inside idxU
		  ii=cntU-1; 
		  // position of the last existing row inside idxU
		  jj=t-1;
		  while (ii>=t && jj>=0) {
	                // currently greatest fill index
			mm=idxU[ii];
			// currently greatest existing index
			tt=idxU[jj];
			// insert 0
			if (mm>tt) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                           *pb=0.0;
#else
                           pb->r=pb->i=0.0;
#endif		      
  		           ii--;
			}
			// we must have mm!=tt by construction! 
			else { // mm<tt, shift existing entry
		           *pb=*pU--;
                           jj--;
		        } // end if-else
			pb--;
	          } // end while
		  while (ii>=t) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                        *pb=0.0;
#else
                        pb->r=pb->i=0.0;
#endif		      
  		        ii--;
			pb--;
	          } // end while
		  while (jj>=0) {
		        *pb--=*pU--;
			jj--;
	          } // end while
	      } // end for i2

	      // copy BUT31 from the previous block and supplement it with zero rows
	      // position of the last entry of BUT31
	      pU=BUT->valE[i_bldu-1]+m*k-1;
	      // futures position of the last entry of the extended BUT31 inside
	      // the enlarged block
	      pb=BUT->valE[i_bldu]+cntU*k-1;
	      // index array column i_bldu-1
	      pi=BUT->rowind[i_bldu-1];
	      // index array column i_bldu
	      pi2=BUT->rowind[i_bldu];
	      // number of rows of BUT21
	      it=BUTfirst[i_bldu-1];
	      for (i2=0; i2<k; i2++,pU-=it) {
	          // position of the last row inside BUT->nblockrow[i_bldu-1]
		  ii=m-1; 
	          // position of the last row inside BUT->nblockrow[i_bldu]
		  jj=t-1;
		  while (ii>=it && jj>=0) {
	                // currently largest fill index
			mm=pi[ii];
			// currently largest existing index
			tt=pi2[jj];
			// insert 0
			if (tt>mm) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                           *pb=0.0;
#else
                           pb->r=pb->i=0.0;
#endif		      
  		           jj--;
			}
			else if (tt<mm) { // copy existing entry
		           *pb=*pU--;
                           ii--;
		        } 
			else { // common index, copy and decrement both counters
		           *pb=*pU--;
                           ii--;
                           jj--;
		        } // end if-else
			pb--;
	          } // end while
		  while (ii>=it) {
		        *pb--=*pU--;
                        ii--;
	          } // end while
		  while (jj>=0) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                        *pb=0.0;
#else
                        pb->r=pb->i=0.0;
#endif		      
  		        jj--;
			pb--;
	          } // end while
	      } // end for i2

	      // 2. update BUT32_new = BUT32+BUT31*L21^T
	      // the columns of BUT32 have to be tailored to those
	      // associated with the nonzero rows of L21
	      // therefore we gather those columns to buff
	      // number of rows in L21
	      ii=BLfirst[i_bldu-1];
	      // starting position of BUT32
	      pU=BUT->valE[i_bldu]+cntU*k;
	      // index array of L21
	      pi=BL->rowind[i_bldu-1];
	      // use buff to hold the copy
	      if (cntU*ii>nbuff) {
		 // increase buffer
		 nbuff=ELBOW_BUFF*(cntU*ii)+BUFF_EXT;
		 buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	      } // end if
	      pb=buff;
	      // leading block diagonal entry used to adjust *pi
	      tt=BUT->colind[i_bldu][0];
	      i2=1;
	      for (jj=0; jj<ii; jj++,pb+=cntU,pi++)
		  // copy column *pi of BUT32 to buff
		  COPY(&cntU, pU+(*pi-tt)*cntU,&i2, pb,&i2);
	      
	      // now we are ready to apply GEMM:
	      // buff =   BUT31*L21^T +   BUT32
	      //      = 1*BUT31*L21^T + 1*buff
	      transa="n"; transb="t";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
              alpha=1.0;   beta=1.0;
#else
              alpha.r=1.0; beta.r=1.0;
	      alpha.i=0.0; beta.i=0.0;
#endif
              if (cntU && ii && k)
		 GEMM(transa,transb, &cntU,&ii,&k, &alpha,
                      BUT->valE[i_bldu],&cntU,
		      BL->valE[i_bldu-1],&l, &beta,
		      buff,&cntU, 1,1);
	      // scatter the columns of buff back to BUT32
	      // index array of L21
	      pi=BL->rowind[i_bldu-1];
	      // buff used as copy
	      pb=buff;
	      for (jj=0; jj<ii; jj++,pb+=cntU,pi++)
		  // copy column *pi of BUT32 to buff
		  COPY(&cntU, pb,&i2, pU+(*pi-tt)*cntU,&i2);
	      
	      // 3. finalization
	      // to finalize the aggregated transposed super-diagonal block,
	      // we extend the column index set as well as the row index set
	      // shift existing indices
	      // old index position last column index
	      pi=BUT->colind[i_bldu]+r-1;
	      // new index position last column index
	      pi2=BUT->colind[i_bldu]+ldb-1;
	      // shift indices
	      for (ii=0; ii<r; ii++)
	          *pi2--=*pi--;
	      // insert column indices from the previous block
	      memcpy(BUT->colind[i_bldu],BUT->colind[i_bldu-1],k*sizeof(integer));
	      BUT->nblockcol[i_bldu]=ldb;
	      // insert new row indices in increasing order
	      // position of the last fill row inside idxU
	      ii=cntU-1; 
	      // position of the last existing row inside idxU
	      jj=t-1;
	      // future position of the last index of the extended index set inside
	      // the enlarged block
	      pi=BUT->rowind[i_bldu]+cntU-1;
	      while (ii>=t && jj>=0) {
	            // currently greatest fill index
	            mm=idxU[ii];
		    // currently greatest existing index
		    tt=idxU[jj];
		    // insert mm
		    if (mm>tt) {
		       *pi=mm;
		       ii--;
		    }
		    else { // shift existing index (here copy it back from idxU)
		       *pi=tt;
                       jj--;
		    } // end if-else
		    pi--;
	      } // end while
	      while (ii>=t) {
	            // currently largest fill index
	            mm=idxU[ii];
		    *pi--=mm;
		    ii--;
	      } // end while
	      while (jj>=0) {
		    // currently largest existing index
		    tt=idxU[jj];
		    *pi--=tt;
		    jj--;
	      } // end while
	      BUT->nblockrow[i_bldu]=cntU;

#ifdef PRINT_INFO0
	      printf("extended BUT{%ld}\n",i_bldu);
	      if (BUT->nblockrow[i_bldu]) {
		 printf("        ");
		 for (jj=0; jj<BUT->nblockcol[i_bldu]; jj++) {
		     printf("%8ld",BUT->colind[i_bldu][jj]);
		 }
		 printf("\n");
		 for (ii=0; ii<BUT->nblockrow[i_bldu]; ii++) {
		     printf("%8ld",BUT->rowind[i_bldu][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     for (jj=0; jj<BUT->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BUT->valE[i_bldu][BUT->nblockrow[i_bldu]*jj+ii]);
		     }
#else
		     for (jj=0; jj<BUT->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BUT->valE[i_bldu][BUT->nblockrow[i_bldu]*jj+ii].r);
		     }
		     printf("\n        ");
		     for (jj=0; jj<BUT->nblockcol[i_bldu]; jj++) {
		         printf("%8.1le",BUT->valE[i_bldu][BUT->nblockrow[i_bldu]*jj+ii].i);
		     }
#endif
		     printf("\n");
		 }
		 printf("\n");
	      }
#endif
	      // --------------------------------------------------------
	      // END compute aggregated transposed upper triangular block
	      // --------------------------------------------------------


	      
	      // ----------------------------------------------------
	      // compute aggregated inverse diagonal block
	      // ----------------------------------------------------	      

	      // [D11_new D12_new]=[D11+D11*U12*D22*L21  -D11*U12*D22]
	      // [D21_new D22    ] [     -D22*L21                D22 ]
	      //                  =[D11+D11*buffC^T*L21  -D11*buffC^T]
	      //                   [     -D22*L21                D22 ]
	      //       ...        =[D11-D12_new*L21      -D12_new]
	      //                   [     -D22*L21            D22 ]
	      	      
	      // 1. shifts & copy

	      // shift current (inverse) diagonal block D22 inside the enlarged block
	      // old position last element
	      pD=BiD->valD[i_bldu]+r*r-1;
	      // recall that in the 1x1 case we do we use 1/D22 rather than D22
	      if (!invert_blocks && r==1) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 *pD=1.0/(*pD);
#else
#ifndef _COMPLEX_SYMMETRIC_
		 pD->r=1.0/pD->r;
#else
		 rval=1.0/FABS(*pD);
		 pD->r= pD->r*rval;
		 pD->i=-pD->i*rval;
#endif
#endif
		 pivots[0]=1;
	      } // end if !invert_blocks and r=1
	      // new position last element
	      pb=BiD->valD[i_bldu]+ldb*ldb-1;
	      // shift and remember different leading dimensions r and ldb
	      for (ii=0; ii<r; ii++,pb-=k)
		  for (jj=0; jj<r; jj++)
		      *pb--=*pD--;
  		      
	      // copy D11 to the leading part of the new enlarged diagonal block
	      // this will become D11_new^0 which has still to be updated
	      // start D11 in the previous block
	      pD=BiD->valD[i_bldu-1];
	      // recall that in the 1x1 case we do we use 1/D11 rather than D11
	      if (!invert_blocks && k==1) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 *pD=1.0/(*pD);
#else
#ifndef _COMPLEX_SYMMETRIC_
		 pD->r=1.0/pD->r;
#else
		 rval=1.0/FABS(*pD);
		 pD->r= pD->r*rval;
		 pD->i=-pD->i*rval;
#endif
#endif
		 pivots[-1]=1;
	      } // end if !invert_blocks and k=1
	      // start D11_new
	      pb=BiD->valD[i_bldu];
	      jj=1;
	      for (ii=0; ii<k; ii++,pD+=k,pb+=ldb)
	          // copy column ii of D11
		  COPY(&k, pD,&jj, pb,&jj);

	      // init new sub-diagonal block D21_new of BiD with 0
	      // beginning of the sub-diagonal block D21_new
	      pD=BiD->valD[i_bldu]+k;
	      for (ii=0; ii<k; ii++,pD+=k)
		  for (jj=0; jj<r; jj++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		      *pD++=0.0;
#else
		      pD->r=pD->i=0.0;
		      pD++;
#endif
                  } // end for jj

	      // init new super-diagonal block D12_new of BiD with 0
	      // beginning of the super-diagonal block
	      pD=BiD->valD[i_bldu]+ldb*k;
	      for (ii=0; ii<r; ii++,pD+=r) 
		  for (jj=0; jj<k; jj++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                      *pD++=0.0;
#else
		      pD->r=pD->i=0.0;
		      pD++;
#endif
                  } // end for jj

              // 2. D12

	      if (invert_blocks) {
		 // compute new super-diagonal block 
		 // D12_new =   -D11*U12*D22
		 //         =   -D11*BUT21^T*D22
		 //         =   -D11*BUT21^T*buff
		 //         =   -D11*buffC^T (1.GEMM)
		 //         = -1*D11*buffC^T + 0*D12_new  using GEMM twice
		 // we need to compute D22^T*BUT21 in order to build the
		 // extended inverse diagonal block
		 // extract the rows of D22 that are associated with the nonzero
		 // rows of BUT21, copy them to buff, then apply GEMM twice
		 // number of rows in BUT21
		 ii=BUTfirst[i_bldu-1];
		 // start of D22 in the extended block
		 pD=BiD->valD[i_bldu]+(ldb+1)*k;
		 // use buff to hold the copy
		 if (ii*r>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*(ii*r)+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 // target buffer
		 pb=buff;
		 // pointer to the array of transposed super-diagonal row indices
		 pi=BUT->rowind[i_bldu-1];
		 // first index of the diagonal block, we need this information in
		 // order to correctly adjust the pointer to the rows of D22
		 tt=BiD->colind[i_bldu][0];
		 for (jj=0; jj<ii; jj++,pb++,pi++)
		     // copy logical row *pi of D22 to buff
		     COPY(&r, pD+(*pi-tt),&ldb, pb,&ii);
	      
/*	      
mexPrintf("buff copy:\n");
pb=buff;	      
for (jj=0; jj<ii; jj++,pb++){
    for (i2=0; i2<r; i2++)
        mexPrintf("%8.1le",pb[ii*i2]);
    mexPrintf("\n");
}
*/
	         // cache the product BUT21^T*D22 in buffC
		 jj=k*r;
		 if (jj>nbuffC) {
		    nbuffC=ELBOW_BUFF*jj+BUFF_EXT;
		    buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
		 } // end if
		 // init buffC with 0s
		 pC=buffC;
		 for (i2=0; i2<jj; i2++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     *pC++=0.0;
#else
		     pC->r=pC->i=0.0;
		     pC++;
#endif
		 } // end for i2
		 // now we are ready to apply GEMM for the first time:
		 // buffC =   BUT21^T*D22
		 //       = 1*BUT21^T*buff + 0*buffC
		 transa="t"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=1.0;   beta=0.0;
#else
		 alpha.r=1.0; beta.r=0.0;
		 alpha.i=0.0; beta.i=0.0;
#endif
		 if (k && r && ii)
		    GEMM(transa,transb, &k,&r,&ii, &alpha,
			 BUT->valE[i_bldu-1],&m, 
			 buff,&ii, &beta,
			 buffC,&k,1,1);
/*	      
mexPrintf("buffC:\n");
pb=buffC;	      
for (jj=0; jj<k; jj++,pb++){
    for (i2=0; i2<r; i2++)
        mexPrintf("%8.1le",pb[k*i2]);
    mexPrintf("\n");
}
*/
		 // now we are ready to apply GEMM for the second time:
		 // D12_new = -1*D11*buffC + 0*D12_new
		 transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=-1.0;   beta=0.0;
#else
		 alpha.r=-1.0; beta.r=0.0;
		 alpha.i= 0.0; beta.i=0.0;
#endif
		 if (k && r)
		    GEMM(transa,transb, &k,&r,&k, &alpha,
			 BiD->valD[i_bldu-1],&k,
			 buffC,&k, &beta,
			 BiD->valD[i_bldu]+ldb*k,&ldb,1,1);
		 
/*	      
mexPrintf("D11:\n");
pb=BiD->valD[i_bldu-1];
for (jj=0; jj<k; jj++,pb++){
    for (i2=0; i2<k; i2++)
        mexPrintf("%8.1le",pb[k*i2]);
    mexPrintf("\n");
}
mexPrintf("-D11 buffC:\n");
pb=BiD->valD[i_bldu]+ldb*k;	      
for (jj=0; jj<k; jj++,pb++){
    for (i2=0; i2<r; i2++)
        mexPrintf("%8.1le",pb[ldb*i2]);
    mexPrintf("\n");
}
*/	      

                 // 3. D11	      

		 // compute new leading diagonal block
		 // D11_new = D11*U12*D22*L21 + D11
		 //         =   -D12_new *L21 +   D11 
		 //         = -1*D12_new *L21 + 1*D11 using GEMM
		 // by post-multiplying with L21 we only need a small part of
		 // D12_new associated with the rows of L21
		 // to do so, copy the required columns of D12_new to buff
		 // number of rows in L21
		 ii=BLfirst[i_bldu-1];
		 // start of D12_new
		 pD=BiD->valD[i_bldu]+ldb*k;
		 // target buffer
		 if (k*ii>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*(k*ii)+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 pb=buff;
		 // pointer to the array of sub-diagonal row indices in L
		 pi=BL->rowind[i_bldu-1];
		 // first index of the diagonal block, we need this information in
		 // order to correctly adjust the pointer to the columns of D12_new
		 tt=BiD->colind[i_bldu][0];
		 i2=1;
		 for (jj=0; jj<ii; jj++,pb+=k,pi++)
		     // copy logical column *pi of D12_new to buff
		     COPY(&k, pD+(*pi-tt)*ldb,&i2, pb,&i2);
		 // now we are ready to apply GEMM to compute D11_new:
		 // D11_new = -1*D12_new*L21 + 1*D11 
		 //         = -1*buff   *L21 + 1*D11 
		 transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=-1.0;   beta=1.0;
#else
		 alpha.r=-1.0; beta.r=1.0;
		 alpha.i= 0.0; beta.i=0.0;
#endif
		 if (k && ii)
		    GEMM(transa,transb, &k,&k,&ii, &alpha,
			 buff,&k,
			 BL->valE[i_bldu-1],&l, &beta,
			 BiD->valD[i_bldu],&ldb, 1,1);
	      

		 // 4. D21

		 // compute new sub-diagonal block 
		 // D21_new =   -D22*L21
		 //         = -1*D22*L21+0*D21_new using GEMM
		 // where L21 refers to the leading block of BL->valE[i_bldu-1]
		 // associated with D22. L21 may not have as many rows as D22,
		 // for this reason we have to extract the associated columns of 
		 // D22 in advance
		 // before applying GEMM, copy the columns of D22 associated with
		 // the rows of L21 to buff
		 // number of rows in L21
		 ii=BLfirst[i_bldu-1];
		 // start of D22 in the enlarged diagonal block
		 pD=BiD->valD[i_bldu]+(ldb+1)*k;
		 // use buff to hold the copy
		 if (r*ii>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*(r*ii)+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 // target buffer
		 pb=buff;
		 // pointer to the array of sub-diagonal row indices
		 pi=BL->rowind[i_bldu-1];
		 // first index of the diagonal block, we need this information in
		 // order to correctly adjust the pointer to the columns of D22
		 tt=BiD->colind[i_bldu][0];
		 i2=1;
		 for (jj=0; jj<ii; jj++,pb+=r,pi++)
		     // copy logical row *pi of D22 to buff
		     COPY(&r, pD+(*pi-tt)*ldb,&i2, pb,&i2);
		 // now we are ready to apply GEMM:
		 // D21_new = -1*D22 *L21+0*D21_new
		 //         = -1*buff*L21+0*D21_new
		 transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=-1.0;   beta=0.0;
#else
		 alpha.r=-1.0; beta.r=0.0;
		 alpha.i= 0.0; beta.i=0.0;
#endif
		 if (r && k && ii)
		    GEMM(transa,transb, &r,&k,&ii, &alpha,
			 buff,&r,
			 BL->valE[i_bldu-1],&l, &beta,
			 BiD->valD[i_bldu]+k,&ldb, 1,1);
	      } // end if invert_blocks
	      else {
		 // compute D21_new=L21*(P11*tril(LD11)) using level-2 BLAS
		 // extract the rows of BL21, copy them to buff, then apply GEMV
		 // use buff to hold the copy
		 ii=k*r;
		 if (ii>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*ii+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 // init target buffer
		 pb=buff;
		 for (jj=0; jj<ii; jj++,pb++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     *pb=0.0;
#else
		     pb->r=pb->i=0.0;
#endif
		 } // end for jj
		 // pointer to the target buffer
		 pb=buff;
		 // pointer to the array of sub-diagonal row indices of BL21
		 pi=BL->rowind[i_bldu-1];
		 // first index of the diagonal block D22, we need this information in
		 // order to correctly adjust the pointer to the rows of D22
		 tt=BiD->colind[i_bldu][0];
		 // start of old BL21
		 pD=BL->valE[i_bldu-1];
		 // number of rows in BL21 (=beginning of B31)
		 ii=BLfirst[i_bldu-1];
		 for (jj=0; jj<ii; jj++,pD++,pi++)
		     // copy logical row *pi of BL21 to buff
		     COPY(&k, pD,&l, pb+(*pi-tt),&r);

		 // interchange columns of buffered L21_new
		 // shift pivots vector back
		 pivots-=k;
		 i2=1;
		 for (mm=0; mm<k; mm++) {
		     jj=pivots[mm];
		     jj--; // switch to C-style
			  
		     // if necessary swap columns jj and mm
		     if (jj!=mm)
		        SWAP(&r,buff+jj*r,&i2,buff+mm*r,&i2);
		 } // end for mm
		 // interchange rows of buffered L21_new
		 // shift pivots vector forward
		 pivots+=k;
		 for (mm=0; mm<r; mm++) {
		     jj=pivots[mm];
		     jj--; // switch to C-style
			  
		     // if necessary swap rows jj and mm
		     if (jj!=mm)
		        SWAP(&k,buff+jj,&r,buff+mm,&r);
		 } // end for mm

		 // multiply by L11
		 transa="n"; 
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=1.0;   beta=0.0;
#else
		 alpha.r=1.0; beta.r=0.0;
		 alpha.i=0.0; beta.i=0.0;
#endif
		 // start of BD11(old)
		 pD=BiD->valD[i_bldu-1];
		 // start of BD21_new
		 pb=BiD->valD[i_bldu]+k;
		 for (mm=0; mm<k; mm++) {
		     // proceed with GEMV in order to multiply with L11
		     ii=k-mm;
		     // recall that L is unit block lower triangular
		     // swap diagonal entry
		     val=pD[k*mm+mm];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     pD[k*mm+mm]=1.0;
#else
		     pD[k*mm+mm].r=1.0;
		     pD[k*mm+mm].i=0.0;
#endif
		     GEMV(transa, &r,&ii, &alpha, buff+mm*r,&r,
			  pD+k*mm+mm,&i2, &beta,  pb+mm*ldb,&i2, 1);
		     // swap diagonal entry back
		     pD[k*mm+mm]=val;
		 } // end for mm

		 // compute D12_new=L11^{-1}P11^T*BUT21^T using TRSV
		 // extract the rows of BUT21, copy them to BD12_new, then apply TRSV
		 // pointer to the target BD12_new
		 pb=BiD->valD[i_bldu]+ldb*k;
		 // pointer to the array of sub-diagonal row indices of BUT21
		 pi=BUT->rowind[i_bldu-1];
		 // first index of the diagonal block D22, we need this information in
		 // order to correctly adjust the pointer to the rows of D22
		 tt=BiD->colind[i_bldu][0];
		 // start of old BUT21
		 pD=BUT->valE[i_bldu-1];
		 // number of rows in BUT21 (=beginning of BUT31)
		 ii=BUTfirst[i_bldu-1];
		 for (jj=0; jj<ii; jj++,pD++,pi++)
		     // copy logical row *pi of BUT21 to the associated column of BD12_new
		     COPY(&k, pD,&m, pb+(*pi-tt)*ldb,&i2);

		 // interchange columns of BD12_new
		 // shift pivots vector back
		 pivots-=k;
		 for (mm=0; mm<k; mm++) {
		     jj=pivots[mm];
		     jj--; // switch to C-style
			  
		     // if necessary swap rows jj and mm
		     if (jj!=mm)
		        SWAP(&r,pb+jj,&ldb,pb+mm,&ldb);
		 } // end for mm

		 // solve with L11
		 pD=BiD->valD[i_bldu-1];
		 for (mm=0; mm<r; mm++)
		     TRSV("l","n","u", &k, pD,&k, pb+mm*ldb,&i2, 1,1,1);
		 		 	 
		 // update extended permutation vector
		 // recall that pivots vector is still shifted back!
		 for (mm=k; mm<ldb; mm++)
		     pivots[mm]+=k;
	      } // end if-else invert_blocks


              // 5. finalization
	      
	      // to finalize the aggregated diagonal block, we extend
	      // the column index set
	      // first shift existing indices
	      // old index position last column index
	      pi=BiD->colind[i_bldu]+r-1;
	      // new index position last column index
	      pi2=BiD->colind[i_bldu]+ldb-1;
	      for (ii=0; ii<r; ii++)
		  *pi2--=*pi--;
	      // now insert indices from the previous block
	      memcpy(BiD->colind[i_bldu],BiD->colind[i_bldu-1],k*sizeof(integer));
	      BiD->nblockcol[i_bldu]=ldb;

#ifdef PRINT_INFO0
	      printf("extended BiD{%ld}\n",i_bldu);
	      printf("        ");
	      for (jj=0; jj<BiD->nblockcol[i_bldu]; jj++) {
		  printf("%8ld",BiD->colind[i_bldu][jj]);
	      }
	      printf("\n");
	      for (ii=0; ii<BiD->nblockcol[i_bldu]; ii++) {
		  printf("%8ld",BiD->colind[i_bldu][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  for (jj=0; jj<BiD->nblockcol[i_bldu]; jj++) {
		      printf("%8.1le",BiD->valD[i_bldu][BiD->nblockcol[i_bldu]*jj+ii]);
		  }
#else
		  for (jj=0; jj<BiD->nblockcol[i_bldu]; jj++) {
		      printf("%8.1le",BiD->valD[i_bldu][BiD->nblockcol[i_bldu]*jj+ii].r);
		  }
		  printf("\n        ");
		  for (jj=0; jj<BiD->nblockcol[i_bldu]; jj++) {
		      printf("%8.1le",BiD->valD[i_bldu][BiD->nblockcol[i_bldu]*jj+ii].i);
		  }
#endif
		  printf("\n");
	      }
	      printf("\n");
#endif
	      
	      // ----------------------------------------------------
	      // END compute aggregated inverse diagonal block
	      // ----------------------------------------------------	      



	      // remove old blocks and adapt BLfirst, BUTfirst
	      free(BiD->colind[i_bldu-1]);BiD->colind[i_bldu-1]=NULL;
	      free(BiD->valD[i_bldu-1]);  BiD->valD[i_bldu-1]=NULL;
	      BiD->nblockcol[i_bldu-1]=0;

	      free(BL->rowind[i_bldu-1]);BL->rowind[i_bldu-1]=NULL;
	      free(BL->colind[i_bldu-1]);BL->colind[i_bldu-1]=NULL;
	      free(BL->valE[i_bldu-1]);  BL->valE[i_bldu-1]=NULL;
	      BL->nblockcol[i_bldu-1]=0;
	      BL->nblockrow[i_bldu-1]=0;
	      BLfirst[i_bldu-1]=0;

	      free(BUT->rowind[i_bldu-1]);BUT->rowind[i_bldu-1]=NULL;
	      free(BUT->colind[i_bldu-1]);BUT->colind[i_bldu-1]=NULL;
	      free(BUT->valE[i_bldu-1]);  BUT->valE[i_bldu-1]=NULL;
	      BUT->nblockcol[i_bldu-1]=0;
	      BUT->nblockrow[i_bldu-1]=0;
	      BUTfirst[i_bldu-1]=0;

	      /*
	      */

	      
	   } // end if jj<=ELBOW_BLOCK*ii...
	   else // reset fill counter
	      cnt_fill=0;


	   
	   // reset idxposL, idxposU
	   for (jj=0; jj<cntL; jj++) {
	       ii=idxL[jj];
	       idxposL[ii]=0;
	   } // end for jj
	   cntL=0;
	   for (jj=0; jj<cntU; jj++) {
	       ii=idxU[jj];
	       idxposU[ii]=0;
	   } // end for jj
	   cntU=0;
#ifdef PRINT_CHECK
	   for (m=0; m<n; m++) {
	     if (idxposL[m]) {
	       printf("12.idxposL[%ld]=%ld\n",m,idxposL[m]);
	       fflush(stdout);
	     } // end if
	   } // end for m	
	   for (m=0; m<n; m++) {
	     if (idxposU[m]) {
	       printf("12.idxposU[%ld]=%ld\n",m,idxposU[m]);
	       fflush(stdout);
	     } // end if
	   } // end for m	
#endif
	   
	} // end if i_bldu and PROGRESSIVE_AGGREGATION
#ifdef _PROFILING_
	time_progressive_aggregation+=omp_get_wtime()-timeBegin;
#endif
	

	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// finally update linked lists
	// recall that BLfirst as well as BUTfirst are already advanced such
	// that their first entry refers to an index beyond block column i
	
	// beginning of the linked list
	it=BLhead[i_bldu];
	// while list is non-empty
	while (it>=0) {
	      // position of the current leading entry in block column "it" of BL
	      m=BLfirst[it];
	      // next block column
	      it_next=BLlist[it];
	      // make sure that we did not reach the end
	      if (m<BL->nblockrow[it]) {
		 // leading row index of block column "it"
		 tt=BL->rowind[it][m];
		 // scalar index tt refers to block ii
		 ii=invblock[tt];
		 // update linked list for BL
		 BLlist[it]=BLhead[ii];
		 BLhead[ii]=it;
	      }
	      // restore next block column
	      it=it_next;
	} // end while
	// insert first nonzero row index of block column i of BL
	BLfirst[i_bldu]=0;
	// block column non-empty
	if (0<BL->nblockrow[i_bldu]) {
	   // get row index
	   t=BL->rowind[i_bldu][0];
	   // associated block
	   tt=invblock[t];
	   // add i to the linked list of block column tt
	   BLlist[i_bldu]=BLhead[tt];
	   BLhead[tt]=i_bldu;
	} // end if
#ifdef PRINT_INFO0
	printf("linked list BL\n");
	for (m=0; m<BL->nblocks; m++)
	    printf("%8ld",BLhead[m]);
	printf("\n");
	for (m=0; m<BL->nblocks; m++)
	    printf("%8ld",BLlist[m]);
	printf("\n");
	for (m=0; m<BL->nblocks; m++)
	    printf("%8ld",BLfirst[m]);
	printf("\n");
	fflush(stdout);
#endif
	
	// beginning of the linked list
	it=BUThead[i_bldu];
	// while list is non-empty
	while (it>=0) {
	      // position of the current leading entry in block column "it" of BUT
	      m=BUTfirst[it];
	      // next block column
	      it_next=BUTlist[it];
	      // make sure that we did not reach the end
	      if (m<BUT->nblockrow[it]) {
		 // leading row index of block column "it"
		 tt=BUT->rowind[it][m];
		 // scalar index tt refers to block ii
		 ii=invblock[tt];
		 // update linked list for BUT
		 BUTlist[it]=BUThead[ii];
		 BUThead[ii]=it;
	      }
	      // restore next block column
	      it=it_next;
	} // end while
	// insert first nonzero row index of block column i of BUT
	BUTfirst[i_bldu]=0;
	// block column non-empty
	if (0<BUT->nblockrow[i_bldu]) {
	   // get first row index
	   t=BUT->rowind[i_bldu][0];
	   // associated block
	   tt=invblock[t];
	   // add i to the linked list of block column tt
	   BUTlist[i_bldu]=BUThead[tt];
	   BUThead[tt]=i_bldu;
	} // end if
#ifdef PRINT_INFO0
	printf("linked list BUT\n");
	for (m=0; m<BUT->nblocks; m++)
	    printf("%8ld",BUThead[m]);
	printf("\n");
	for (m=0; m<BUT->nblocks; m++)
	    printf("%8ld",BUTlist[m]);
	printf("\n");
	for (m=0; m<BUT->nblocks; m++)
	    printf("%8ld",BUTfirst[m]);
	printf("\n");
	fflush(stdout);
#endif
#ifdef _PROFILING_
	time_list_update+=omp_get_wtime()-timeBegin;
#endif
	
	// advance permutation to the next block
	if (!invert_blocks)
	   pivots+=BiD->nblockcol[i_bldu];

	// advance to next block column
	i++;
	i_bldu++;
  } // end while j<n
  // **************************************************************************
  // *****                     END main loop                              *****
  // **************************************************************************
  if (!invert_blocks)
     pivots-=n;


#ifdef _PROFILING_
  timeBegin=omp_get_wtime();
#endif
  // post processing

  // shift blocks in order to erase old empty blocks
  i_bldu=0;
  for (i=0; i<BiD->nblocks; i++) {
      // block size of the current block
      bi=BiD->nblockcol[i];

      // if the block was non-empty
      if (bi!=0) {
	 // shift block structures
	 BiD->nblockcol[i_bldu]=BiD->nblockcol[i];
	 BiD->colind[i_bldu]   =BiD->colind[i];
	 BiD->valD[i_bldu]     =BiD->valD[i];

	 BL->nblockcol[i_bldu]=BL->nblockcol[i];
	 BL->nblockrow[i_bldu]=BL->nblockrow[i];
	 BL->colind[i_bldu]   =BL->colind[i];
	 BL->rowind[i_bldu]   =BL->rowind[i];
	 BL->valE[i_bldu]     =BL->valE[i];
	 
	 BUT->nblockcol[i_bldu]=BUT->nblockcol[i];
	 BUT->nblockrow[i_bldu]=BUT->nblockrow[i];
	 BUT->colind[i_bldu]   =BUT->colind[i];
	 BUT->rowind[i_bldu]   =BUT->rowind[i];
	 BUT->valE[i_bldu]     =BUT->valE[i];

	 i_bldu++;
      } // end if
  } // end for i
  BiD->nblocks=BL->nblocks=BUT->nblocks=i_bldu;
  
  // re-arrange sparse block data structures
  BiD->nblockcol=(integer *) realloc(BiD->nblockcol,(size_t)BiD->nblocks*sizeof(integer));
  BiD->colind   =(integer **)realloc(BiD->colind,   (size_t)BiD->nblocks*sizeof(integer *));
  BiD->valD     =(FLOAT **)  realloc(BiD->valD,     (size_t)BiD->nblocks*sizeof(FLOAT *));

  BL->nblockcol=(integer *) realloc(BL->nblockcol,(size_t)BL->nblocks*sizeof(integer));
  BL->nblockrow=(integer *) realloc(BL->nblockrow,(size_t)BL->nblocks*sizeof(integer));
  BL->colind   =(integer **)realloc(BL->colind,   (size_t)BL->nblocks*sizeof(integer *));
  BL->rowind   =(integer **)realloc(BL->rowind,   (size_t)BL->nblocks*sizeof(integer *));
  BL->valE     =(FLOAT **)  realloc(BL->valE,     (size_t)BL->nblocks*sizeof(FLOAT *));

  BUT->nblockcol=(integer *) realloc(BUT->nblockcol,(size_t)BUT->nblocks*sizeof(integer));
  BUT->nblockrow=(integer *) realloc(BUT->nblockrow,(size_t)BUT->nblocks*sizeof(integer));
  BUT->colind   =(integer **)realloc(BUT->colind,   (size_t)BUT->nblocks*sizeof(integer *));
  BUT->rowind   =(integer **)realloc(BUT->rowind,   (size_t)BUT->nblocks*sizeof(integer *));
  BUT->valE     =(FLOAT **)  realloc(BUT->valE,     (size_t)BUT->nblocks*sizeof(FLOAT *));
  
  
  
  // recompute block diagonal scaling for U and now store it there
  for (i=0; i<BiD->nblocks; i++) {
      // multiply (transposed) strict upper triangular block by the inverse
      // diagonal block buffer it and copy it back to BUT
      pU=BUT->valE[i];
      bi=BUT->nblockcol[i];
      l =BUT->nblockrow[i];
      // scalar case, skip LAPACK
      if (bi==1) {
	 // inverse diagonal entry
	 val=BiD->valD[i][0];
	 for (ii=0; ii<l; ii++,pU++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     *pU*=val;
#else
	     rval =pU->r*val.r-pU->i*val.i;
	     pU->i=pU->r*val.i+pU->i*val.r;
	     pU->r=rval;
#endif
	 } // end for ii
      }
      else {
	 k=l*bi;    // total size
	 if (k>nbuff) {
	    // increase buffer
	    nbuff=ELBOW_BUFF*k+BUFF_EXT;
	    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	 } // end if
	 if (invert_blocks) {
	    // use level-3-BLAS 
	    // init with 0
	    pC=buff;
	    for (ii=0; ii<k; ii++,pC++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	        *pC=0.0;
#else
		pC->r=pC->i=0.0;
#endif
	    } // end for ii
	    pC=buff;
	    // C = 1*A*B^T+0*C
	    transa="n"; transb="t";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    alpha=1.0;
	    beta =0.0;
#else
	    alpha.r=1.0; alpha.i=0.0;
	    beta.r =0.0; beta.i =0.0;
#endif
	    if (l && bi)
	       GEMM(transa,transb, &l,&bi,&bi, &alpha,
		    pU,&l,
		    BiD->valD[i],&bi, &beta, pC,&l, 1,1);
	    // copy values back to block column i of BUT
	    // pU <- pC
	    memcpy(pU,pC,(size_t)k*sizeof(FLOAT));
	 } // end if invert_blocks
	 else { // LU case
	    // forward/backward solve using the LU factorization of the diagonal block
	    // copy subdiagonal part of buff to BUT since solve works in-place with BUT^T
	    // buff = BUT(:,i)^T
	    pb=buff;
	    i2=1;
	    for (ii=0; ii<l; ii++,pb+=bi,pU++) 
	        // pU -> pb
	        COPY(&bi,pU,&l,pb,&i2);
	    // forward/backward solves LU * BUT^T  <- BUT^T
	    //                    <=>  LU * buffC <- buffC
	    GETRS("n", &bi, &l, BiD->valD[i],&bi, pivots, buff,&bi, &ierr,1);
	    if (ierr<0) {
	       printf("LAPACK's GETRS: %ld-th argument had an illegal value\n", -ierr);
	       return (ierr);
	    }
	    // BUT(:,i)=buff^T
	    pb=buff;
	    pU=BUT->valE[i];
	    i2=1;
	    for (ii=0; ii<l; ii++,pb+=bi,pU++) 
	        // pb -> pU
	        COPY(&bi,pb,&i2,pU,&l);	      
	   } // end if-else invert_blocks
      } // end if-else bi==1
      
      // advance permutation to the next block
      if (!invert_blocks)
	 pivots+=bi;
  } // end for i
  
#ifdef PRINT_INFO0
  if (!invert_blocks)
     pivots-=n;
  
  printf("Final factorization\n");
  printf("BiD, BL, BUT\n");fflush(stdout);
  for (i=0; i<BiD->nblocks; i++) {
      printf("block %3ld of BL, BiD, BUT\n",i);fflush(stdout);
      if (invert_blocks)
	 printf("inverse diagonal block BiD[%3ld]\n        ",i);
      else
	 printf("factorized diagonal block BiD[%3ld]\n        ",i);
      for (m=0; m<BiD->nblockcol[i]; m++) 
	  printf("%8ld",BiD->colind[i][m]);
      printf("\n");fflush(stdout);
      for (k=0; k<BiD->nblockcol[i]; k++) {
	  printf("%8ld",BiD->colind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%8.1le",BiD->valD[i][m*BiD->nblockcol[i]+k]);
	  printf("\n");fflush(stdout);
#else	    
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%8.1le",BiD->valD[i][m*BiD->nblockcol[i]+k].r);
	  printf("\n        ");
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%8.1le",BiD->valD[i][m*BiD->nblockcol[i]+k].i);
	  printf("\n");fflush(stdout);
#endif
      }
      printf("\n");fflush(stdout);
      if (!invert_blocks) {
	 printf("pivots: ");
	 for (m=0; m<BiD->nblockcol[i]; m++)
	     printf("%8ld",pivots[m]);
	 printf("\n");fflush(stdout);
	 pivots+=BiD->nblockcol[i];
      }

      printf("sub-diagonal block BL[%3ld]\n        ",i);
      if (BL->nblockrow[i]>0) {
	 for (m=0; m<BL->nblockcol[i]; m++) 
	     printf("%8ld",BL->colind[i][m]);
	 printf("\n");fflush(stdout);
	 for (k=0; k<BL->nblockrow[i]; k++) {
	     printf("%8ld",BL->rowind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%8.1le",BL->valE[i][m*BL->nblockrow[i]+k]);
#else
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%8.1le",BL->valE[i][m*BL->nblockrow[i]+k].r);
	     printf("\n        ");
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%8.1le",BL->valE[i][m*BL->nblockrow[i]+k].i);
#endif
	     printf("\n");fflush(stdout);
	 } // end for k
      }
      printf("\n");fflush(stdout);

      printf("transposed super-diagonal block BUT[%3ld]\n        ",i);
      if (BUT->nblockrow[i]>0) {
	 for (m=0; m<BUT->nblockcol[i]; m++) 
	     printf("%8ld",BUT->colind[i][m]);
	 printf("\n");fflush(stdout);
	 for (k=0; k<BUT->nblockrow[i]; k++) {
	     printf("%8ld",BUT->rowind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     for (m=0; m<BUT->nblockcol[i]; m++) 
	         printf("%8.1le",BUT->valE[i][m*BUT->nblockrow[i]+k]);
#else
	     for (m=0; m<BUT->nblockcol[i]; m++) 
	         printf("%8.1le",BUT->valE[i][m*BUT->nblockrow[i]+k].r);
	     printf("\n        ");
	     for (m=0; m<BUT->nblockcol[i]; m++) 
	         printf("%8.1le",BUT->valE[i][m*BUT->nblockrow[i]+k].i);
#endif
	     printf("\n");fflush(stdout);
	 } // end for k
      }
      printf("\n");fflush(stdout);
  } // end for i
#endif

  
  // sort original matrix back such that row indices and their
  // values are taken in increasing order
  for (j=0; j<n; j++)
      QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));


  // release memory
  idxL--;free(idxL);
  free(idxposL);
  idxU--;free(idxU);
  free(idxposU);
  free(Ahead);
  free(Alist);
  free(Afirst);
  free(BLhead);
  free(BLlist);
  free(BLfirst);
  free(BUThead);
  free(BUTlist);
  free(BUTfirst);
  
  free(AvalD);
  free(AvalL);
  free(AvalU);
  free(Arowind);
  free(Acolind);
  free(buffC);
  free(buff);

  // no a priori block partitioning given, use ILU(1,droptol)
  if (flag_ilu1t) {
     free(blocksize);
     // restore user's values
     blocksize=myblocksize;
     nblocks  =mynblocks;
  } // end if flag_ilu1t

  if (perturbation)
     buffDGLBLCK=FREE(buffDGLBLCK);

  // finally include diagonal scaling to the determinant
  if (SL!=NULL)
     for (i=0; i<n; i++)
         determinant->r-=log(SL[i]);
  if (SR!=NULL)
     for (i=0; i<n; i++)
         determinant->r-=log(SR[i]);

  // inverse permutation w.r.t. p
  for (i=0; i<n; i++) {
      invblock[i]=p[i];
      idxbuff[invblock[i]]=i;
  } // end for i
  // compute sign of p
  flag=0;
  for (i=0; i<n; i++) {
      k=invblock[i];
      if (k!=i) {
	 // at which position is i located inside p?
	 j=idxbuff[i];
	 // swap components i and j of p
	 invblock[i]=invblock[j]; invblock[j]=k;
	 // swap components i and k of p^{-1}
	 idxbuff[i]=idxbuff[k]; idxbuff[k]=j;
	 // toggle flag to compute the sign of the permutation
	 flag=!flag;
      } // end if
  } // end for i
  // sign change by permutation
  if (flag)
     determinant->i+=M_PI;

  // inverse permutation w.r.t. invq
  for (i=0; i<n; i++) {
      invblock[i]=invq[i];
      idxbuff[invblock[i]]=i;
  } // end for i
  // compute sign of invq
  flag=0;
  for (i=0; i<n; i++) {
      k=invblock[i];
      if (k!=i) {
	 // at which position is i located inside invq?
	 j=idxbuff[i];
	 // swap components i and j of invq
	 invblock[i]=invblock[j]; invblock[j]=k;
	 // swap components i and k of invq^{-1}
	 idxbuff[i]=idxbuff[k]; idxbuff[k]=j;
	 // toggle flag to compute the sign of the permutation
	 flag=!flag;
      } // end if
  } // end for i
  // sign change by permutation
  if (flag)
     determinant->i+=M_PI;
  
  free(invblock);
  free(idxbuff);

  // make sure to consider the main branch
  while (determinant->i<=-M_PI)
        determinant->i+=2*M_PI;
  while (determinant->i>M_PI)
        determinant->i-=2*M_PI;

  BL->nnz=0;
  for (i=0; i<BL->nblocks-1; i++) 
      BL->nnz+=BL->nblockrow[i]*BL->nblockcol[i];
  BiD->nnz=0;
  for (i=0; i<BiD->nblocks; i++) 
      BiD->nnz+=BiD->nblockcol[i]*BiD->nblockcol[i];
  BUT->nnz=0;
  for (i=0; i<BUT->nblocks-1; i++) 
      BUT->nnz+=BUT->nblockrow[i]*BUT->nblockcol[i];
  
#ifdef _PROFILING_
  time_post_processing+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
  time_bilu=omp_get_wtime()-time_bilu;
  printf("profiling summary\n");
  printf("0) initialization matrix and lists            %8.1le\n",time_init_matrix);
  if (ilu1t==BLOCK_ILU1T)
     printf("1) ILU(1,%8.1le)                            %8.1le\n",droptol,time_ilu1t);
  else if (ilu1t==BLOCK_ILUPT)
     printf("1) ILU(%2ld,%8.1le)                           %8.1le\n",level_of_fill,droptol,time_ilu1t);
  else if (ilu1t==BLOCK_SUPERNODES)
     printf("1) SUPERNODES                                 %8.1le\n",time_ilu1t);
  else
     printf("1) ---                                        %8.1le\n",time_ilu1t);
  printf("2) Extraction of blocks from A                %8.1le\n",time_extract_block);
  printf("3) LU update first pass                       %8.1le\n",time_LU_update_pass_1);
  printf("   Scalar LU update                           %8.1le\n",time_scalar_update);
  printf("   Rank-1 LU update                           %8.1le\n",time_rank_1_update);
  printf("   Small size block LU update                 %8.1le\n",time_small_block_update);
  printf("   Level-3 BLAS LU gather                     %8.1le\n",time_level_3_blas_gather);
  printf("   Level-3 BLAS LU update                     %8.1le\n",MAX(0.0,time_level_3_blas_update-time_level_3_blas_gather-time_level_3_blas_scatter));
  printf("   Level-3 BLAS LU scatter                    %8.1le\n",time_level_3_blas_scatter);
  printf("   Linked list update                         %8.1le\n",time_list_update);
  printf("   Remainder LU update pass 2                 %8.1le\n",time_LU_update_pass_2);
  printf("4) Diagonal block factorization and inversion %8.1le\n",time_diagonal_block);
  printf("5) Off-diagonal block re-scaling and dropping %8.1le\n",time_off_diagonal_block);
  printf("6) Progressive aggregation                    %8.1le\n",time_progressive_aggregation);
  printf("7) Post processing                            %8.1le\n",time_post_processing);
  printf("Total BILU time %8.1le\n\n",time_bilu);

  fflush(stdout);
#endif
  
  return (ierr);
} // end bilu

 

 
