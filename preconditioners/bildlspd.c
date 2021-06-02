/* $Id: bildlspd.c 7315 2021-05-28 21:00:20Z bolle $

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	September 10, 2017. JANUS Block ILU R1.0.  

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


#define PROGRESSIVE_AGGREGATION 1
// #define PROGRESSIVE_AGGREGATION 0

// #define PRINT_CHECK
// #define PRINT_INFO
// #define printf mexPrintf
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

integer BILDLSPD(SPARSEMATRIX *A,
		 SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD,
		 SPARSEBLOCKMATRIX *BUT,
		 REALS *SL, REALS *SR, integer *p, integer *invq,
		 doublecomplex *determinant,
		 integer *blocksize, integer nblocks,
		 REALS droptol, integer ildl1t, integer pa,
		 integer perturbation, integer invert_blocks,
		 integer level_of_fill)
{  
  /*
    compute block incomplete LDL^T factorization with relative drop tolerance

    This routine approximately computes P^T S_L A S_L P ~  L D L^T,
    where P is a given permutation matrix, S_L is a real 
    diagnal scaling matrix, 
    L is lower unit block triangular matrix and
    D is a block diagonal matrix
    matrix.

    Input
    -----
    A           assumed to be a nonsingular sparse SPD matrix such that
                only half of the matrix is stored (e.g. lower triangular part)
    p,invq      permutation p and inverse permutation invq associated with the a priori
                permutation matrix P
    SL          real diagonal scaling matrices
    SR          not referenced
    blocksize   array with the initial block sizes of each diagonal block
                of the scaled and reordered system.
    nblocks     length of array blocksize		
    perturbation if different from zero, then the diagonal part is modified
                 for positive perturbation the Ajiz/Jennings strategy is used
                 for negative perturbation the "minus" Ajiz/Jennings strategy is employed

    Output
    ------
    BL          lower block triangular matrix L with unit diagonal
    BiD         block diagonal matrix D^{-1} (LAPACK-factorized and inverted)
    BUT         not referenced
  */

  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,mm,tt,it,it_next,jt,jjt,jt_next,i2, // counters
          *invblock=NULL,// mapping scalar index j -> block i of A
          nrAL,ncAL,nrncAL,// estimates for the blocks of A
          *pi,*pi2,    // index array pointers
          nbuff,nbuffC,nbuffD,// logical size of buff, buffC, buffD
          lda,ldb,ldc, // leading dimension used in combination with level-3-BLAS
          flag,        // flag indicating when a block column is exceeded
          bi,          // current block size
          nrowind,     // length of array Arowind
          i_bldl,      // index counter for the block ILDL
          cnt_fill,    // counter for additional fill when merging blocks
          *rowind,ncol,*colind,nrow,// short cuts
    
          *idxL=NULL,*idxU=NULL,// temporary index list of row indices
          *idxposL=NULL,*idxposU=NULL,// associated inverse mapping storing the position
          cntL,        // length of idxL
          cntL_old,    // old length of idxL when needed

          cntU,        // length of idxU
          cntU_old,    // old length of idxU when needed

          *idxbuff=NULL,// auxiliary index buffer
          skip_gthr_sctr,// flag when gather/scatter operation is not needed

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
    
          *Arowind=NULL,// array with the nonzero row indices of the current
                       // sub-diagonal block of A

          flag_ildl1t=0,// flag for using ILDL(1,droptol)
          mynblocks,   // copies of nblocks and blocksize, when ILDL(1,droptol)
          *myblocksize,// is called, copied back on return
    
          ierr=0;      // return value bilu

  
  FLOAT   **AvalD=NULL,// buffer for the diagonal blocks of A
          *AvalL=NULL, // buffer for the current sub-diagonal block of A
          *AvalU=NULL, // buffer for the current transposed super-diagonal block of A
          *pL,*pD,*pU,*pb,*pC, // pointers to lower/diagonal/upper/buff/buffC
          val,valD,    // temporary scalar numerical value
          *buff=NULL,  // buffer for the current block column/row of BL
          *buffC=NULL, // buffer for level-3-BLAS update
          *buffD=NULL, // buffer for diagonal block x upper triangular block
          *buffDgl=NULL,//keep track on the diagonal entries
          ajj,         // auxiliary value for computing the determinant
          alpha, beta; // parameters for level-3-BLAS

  doublecomplex locdet;
  
  REALS   rval,srval,sb,e,ab,d;  // temporary scalar real numerical value

  char    *transa, *transb; // strings for level-3-BLAS

#ifdef _PROFILING_
  double timeBegin,
         time_bilu=0.0,
         time_init_matrix=0.0,
         time_ilu1t=0.0,
         time_extract_block=0.0,
         time_LDL_update_pass_1=0.0,
         time_LDL_update_pass_2=0.0,
         time_scalar_update=0.0,
         time_rank_1_update=0.0,
         time_small_block_update=0.0,
         time_level_3_blas_update=0.0,
         time_list_update=0.0,
         time_diagonal_block=0.0,
         time_off_diagonal_block=0.0,
         time_progressive_aggregation=0.0,
         time_post_processing=0.0;

  time_bilu=omp_get_wtime();
#endif
  
#ifdef PRINT_INFO
  printf("start BICHOL\n");fflush(stdout);
#endif
  
  determinant->r=0.0;
  determinant->i=0.0;

  // direct solver case
  if (droptol<=0.0) {
     droptol=0.0;
     perturbation=0;
  } // end if
  
  // create auxiliary index buffer
  idxbuff=(integer *)malloc((size_t)n*sizeof(integer));
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc((size_t)(n+1)*sizeof(integer));
  idxL[0]=0; idxL++;
  // init with zeros
  idxposL=(integer *)calloc(n,sizeof(integer));
  // auxiliary index buffer
  idxU   =(integer *)malloc((size_t)(n+1)*sizeof(integer));
  idxU[0]=0; idxU++;
  // init with zeros
  idxposU=(integer *)calloc(n,sizeof(integer));

  // buffer for the current block column of BL
  nbuff =n;
  buff  =(FLOAT *)malloc((size_t)nbuff *sizeof(FLOAT));
  // buffer for level-3-BLAS
  nbuffC=n;
  buffC =(FLOAT *)malloc((size_t)nbuffC*sizeof(FLOAT));
  // buffer for diagonal block x upper triangular part
  nbuffD=n;
  buffD =(FLOAT *)malloc((size_t)nbuffD*sizeof(FLOAT));
  // buffer for the diagonal entries
  if (perturbation)
     buffDgl =(FLOAT *)malloc((size_t)n*sizeof(FLOAT));


  
#ifdef PRINT_CHECK
  for (j=0; j<n; j++) {
      if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n)
	 printf("BICHOL ERROR, wrong p/invq, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j], invq[p[j]]);
  }  
  if (SL!=NULL) {
     for (j=0; j<n; j++) {
       if (fabs(SL[j])==0.0)
	  printf("BICHOL ERROR, zero scaling, step %ld\n", j);
     }  
  }
  for (j=0; j<n; j++) {
      for (k=0; k<A->ncol[j]; k++) {
          i=A->rowind[j][k];
	  if (i<0 || i>=n)
	     printf("BICHOL ERROR, illegal matrix index (%ld,%ld)\n", i,j);
     }  
  }
#endif
  // sort columns of A(p,p) such that for every column
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
      for (k=0; k<ncol; k++)
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
#ifdef PRINT_INFO
  printf("system reordered\n");fflush(stdout);
#endif

  // extract diagonal part of P^T*SL*A*SL*P
  if (perturbation) {
     for (j=0; j<n; j++) {
         jj=p[j]; // short cut
	 k=A->ncol[jj]; // nz of column A(:,p[j])
	 rowind=A->rowind[jj];
	 // j-th diagonal entry of P^T*SL*A*SL*P initialized with 0
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	 buffDgl[j]=0.0;
#else
	 buffDgl[j].r=buffDgl[j].i=0.0;
#endif
         for (tt=0; tt<k; tt++) {
	     // index ii of A(ii,jj)
	     ii=rowind[tt];
	     // index p[i]=ii, i.e., A(ii,jj)=A(p[i],p[j])
	     // i=invq[ii];
	     if (ii==jj) {
	        val=A->val[jj][tt];
		if (SL!=NULL) {
		   rval=SL[ii]*SL[jj];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		   val*=rval;
#else
		   val.r*=rval;
		   val.i=0.0; // Hermitian matrix has real diagonal entries
#endif
		} // end if SL!=0
		buffDgl[j]=val;
		// skip column jj
	        tt=k;
	     } // ii==jj
	 } // end for tt
     } // end for j
  } // end if
#ifdef PRINT_INFO
  printf("diagonal part extracted\n");fflush(stdout);
#endif


  /*
  for (j=0; j<n; j++) {
      printf("%12.4le",buffDgl[j]);
  } // end for j
  printf("\n");fflush(stdout);
  */
  
  
  // set up arrays for A(p,p)^T (A is stored by by columns)
  // implicitly these structures work with A(p,p)
  // head of the linked list for the rows of A(p,p)
  Ahead =(integer *)malloc((size_t)n*sizeof(integer));
  // linked list for the leading rows of A(p,p) up to step j
  Alist =(integer *)malloc((size_t)n*sizeof(integer));
  // position of the first entry at step j
  Afirst=(integer *)malloc((size_t)n*sizeof(integer));
  // clear head of the linked list (empty) 
  for (j=0; j<n; j++)
      Ahead[j]=-1;
  // init structures
  for (j=0; j<n; j++) {
      jj=p[j]; // short cut
      // we have to make sure that we are still inside column p[j]
      if (0<A->ncol[jj]) {
	 // first row index k of A(p[k],p[j])
	 k=invq[A->rowind[jj][0]];
	 // position of the current first nonzero entry in column j
	 Afirst[j]=0;
	 // add new entry to the head of the list
         Alist[j]=Ahead[k];
	 Ahead[k]=j;
      } // end if
  } // end for j
#ifdef PRINT_INFO
  printf("linked list set up\n");fflush(stdout);
#endif
#ifdef _PROFILING_
  time_init_matrix=omp_get_wtime()-timeBegin;
#endif
  /*
  for (i=0; i<n; i++)
    printf("%12ld",Ahead[i]);
  printf("\n");
  for (i=0; i<n; i++)
    printf("%12ld",Alist[i]);
  printf("\n");
    for (i=0; i<n; i++)
    printf("%12ld",Afirst[i]);
  printf("\n");
  */

  // no a priori block partitioning given, use ILDL(level_of_fill,droptol)
#ifdef _PROFILING_
  timeBegin = omp_get_wtime();
#endif
  if (nblocks==0 || blocksize==NULL) {
     flag_ildl1t=-1;
     myblocksize=blocksize;
     mynblocks=nblocks;
     blocksize=(integer *)malloc((size_t)n*sizeof(integer));
     if (ildl1t) {
        if (ildl1t==BLOCK_ILU1T) {
#ifdef PRINT_INFO
	   printf("call ILDL1T\n");fflush(stdout);
#endif
	   ierr=ILDL1T_BLOCKS(A, SL,SL, p,invq, blocksize,&nblocks, droptol);
	}
	else if (ildl1t==BLOCK_ILUPT) {
#ifdef PRINT_INFO
	   printf("call ILDLPT\n");fflush(stdout);
#endif
	   ierr=ILDLPT_BLOCKS(A, SL,SR, p,invq, blocksize,&nblocks, droptol,level_of_fill);
	}
	else if (ildl1t==BLOCK_SUPERNODES) {
#ifdef PRINT_INFO
	   printf("call Supernode partitioning\n");fflush(stdout);
#endif
	   ierr=SUPERNODES(A, SL,SL, p,invq, blocksize,&nblocks, droptol);
	}
	else if (ildl1t==BLOCK_APP_SUPERNODES) {
#ifdef PRINT_INFO
	   printf("call approximate supernode partitioning\n");fflush(stdout);
#endif
	   ierr=SPDAPPSUPERNODES(A, SL,SL, p,invq, blocksize,&nblocks, droptol);
	}
     }
     else {
        nblocks=n; for (i=0;i<n;i++)blocksize[i]=1;ierr=0;
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
	buffD    =FREE(buffD    );
        Ahead    =FREE(Ahead    );
        Alist    =FREE(Alist    );
        Afirst   =FREE(Afirst   );
        blocksize=FREE(blocksize);
	
	if (perturbation)
	   buffDgl=FREE(buffDgl);
	
        return (ierr);
     }
  } // end if nblocks=0 or blocksize=0
#ifdef PRINT_CHECK
  for (j=0; j<n; j++) {
      if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n)
	 printf("BICHOL ERROR, wrong p/invq after ILDL1T, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j], invq[p[j]]);
  }  
  if (nblocks>0 && blocksize!=NULL) {
     k=0;
     for (j=0; j<nblocks; j++) {
         k+=blocksize[j];
         if (blocksize[j]<=0 || blocksize[j]>n)
	    printf("BICHOL ERROR, wrong block size, step %ld, blocksize[%ld]=%ld\n", j,j,blocksize[j]);
         if (k>n)
	    printf("BICHOL ERROR, cummulative block size exceeds system size, step %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
     }  
     if (k!=n)
        printf("BICHOL ERROR, cummulative block does not match system size, nblocks %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
     /*
     for (j=0; j<nblocks; j++)
        printf("%4ld",blocksize[j]);
     printf("\n");
     */
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
  
  // set up sparse block data structures for BL,BiD
  BL->nc=BL->nr=n;
  BL->nblocks =nblocks;
  BL->nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BL->nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BL->colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BL->rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BL->valD     =NULL;
  BL->valE     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));

  BiD->nc=BiD->nr=n;
  BiD->nblocks =nblocks;
  BiD->nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
  BiD->nblockrow=NULL;
  BiD->colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
  BiD->rowind   =NULL;
  BiD->valD     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));
  BiD->valE     =NULL;
  for (i=0; i<nblocks; i++) {
      BL->nblockrow[i]=0;
      BL->nblockcol[i]=0;
      BL->colind[i]=BL->rowind[i]=BiD->colind[i]=NULL;
      BL->valE[i]=BiD->valD[i]=NULL;
  }

  // set up arrays for the rows BL^T (BL is stored by block columns)
  // these are required to access BL by rows during the computation
  // of the current block column of BL
  // head of the linked list for the rows of BL
  BLhead =(integer *)malloc((size_t)nblocks*sizeof(integer));
  // linked list for the leading rows of BL up to step j
  BLlist =(integer *)calloc(nblocks,sizeof(integer));
  // position of the first entry at step j
  BLfirst=(integer *)calloc(nblocks,sizeof(integer));
  // now clear head of the linked list (empty) 
  for (j=0; j<nblocks; j++)
      BLhead[j]=-1;

  

  
  // buffer for the current diagonal block, block column/block row of A(p,p)
  nrAL=0;
  ncAL=0;
  nrncAL=0;
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
		 i2=invblock[t]; // i2<=t<r, we may use idxposL as well
		 // new nonzero row index found in this block column
		 if (!idxposU[i2]) {
		    idxU[cntU]=i2;
		    idxposU[i2]=++cntU;
		    (BL->nblockrow[i2])++;
		 } // end if
	      } // end if t<r
	  } // end for m
	  // (partially) reset idxposU for the indices above the diagonal block
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
  } // end for i
  for (i=0; i<nblocks; i++) {
      // printf("%ld\n",BL->nblockrow[i]);fflush(stdout);
      bi=blocksize[i];
      // maximum number of sub-diagonal indices
      nrAL=MAX(nrAL,BL->nblockrow[i]);
      // maximum block size sub-diagonal block
      nrncAL=MAX(nrncAL,bi*BL->nblockrow[i]);
  } // end for i
  // printf("%ld,%ld\n",nrncAL,nrAL);fflush(stdout);
  // finally we have a relatively precise estimate how big the buffers should be
  AvalL     =(FLOAT *)  malloc((size_t)nrncAL *sizeof(FLOAT));
  Arowind   =(integer *)malloc((size_t)nrAL   *sizeof(integer));
  // use array of pointers for the diagonal blocks since we need Dii and
  // Dii^{-1} simultaneously
  AvalD     =(FLOAT  **)malloc((size_t)nblocks*sizeof(FLOAT *));
  for (i=0; i<nblocks; i++)
      AvalD[i]=NULL;
#ifdef _PROFILING_
  time_init_matrix+=omp_get_wtime()-timeBegin;
#endif
  
  
  // **************************************************************************
  // *****                       main loop                                *****
  // **************************************************************************
  i=0; // index for the blocks of A
  j=0; // index for the scalar columns/rows of A
  i_bldl=0; // index for the block of the ILDL, 0<=i_bldl<=i
  cntL=cntU=0; // counter for temporary index lists
  cnt_fill=0; // counter for additional fill-in when merging blocks
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
  while (j<n) {

        // store current block size for easy access
	bi=blocksize[i];

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// temporarily store each diagonal block since we still need
	// it for the Schur complement downdates
	AvalD[i_bldl]=(FLOAT *)malloc((size_t)bi*bi*sizeof(FLOAT));
        // clear buffers
        pD=AvalD[i_bldl];
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
	// for (k=0; k<nrAL; k++) 
	//     Arowind[k]=0;
#ifdef _PROFILING_
	time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	        

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// ----------------------------------------------
	// ----------------------------------------------
	// ------- extract COLUMNS of A(p,p) of block i ------
	// pointer to the current column inside diagonal block i
        pD=AvalD[i_bldl];
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
	BL->colind[i_bldl] =(integer *)malloc((size_t)bi*sizeof(integer));
	BiD->colind[i_bldl]=(integer *)malloc((size_t)bi*sizeof(integer));
	BL->nblockcol[i_bldl]=BiD->nblockcol[i_bldl]=bi;
	// number of extracted sub-blockdiagonal rows
	nrowind=0;
	// loop through the columns r:s-1 inside block column i
	for (k=0; k<bi; k++,pD+=bi,pL+=BL->nblockrow[i_bldl]) {
	    // trivial column indices 
	    BL->colind[i_bldl][k]=BiD->colind[i_bldl][k]=j;

	    jj=p[j]; // shortcut
	    // check COLUMN j of A(p,p)
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
		   rval=SL[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		   val*=rval;
#else
		   val.r*=rval;
		   val.i*=rval;
#endif
		} // end if SL!=0
		
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
		      // nrAL+=10;
		      // Arowind=(integer *)realloc(Arowind,(size_t)nrAL*sizeof(integer));
		      // }
		      // bookmark index for later purposes
		      // Arowind[nrowind++]=t;
		      nrowind++;
		      // nrowind and cntL are identical; similarly, idxL and Arowind coincide
		   } // end if
		   // whoops, our estimate about the number of nonzeros in the
		   // current block row was too small (no idea how this could
		   // have happened)
		   if (nrowind>BL->nblockrow[i]) {
		      // printf("whoops, step %ld, A, realloc\n",i);
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
			  } // end for mm
		      } // end for jj
		      // set additional element in column 0 to 0.0
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
	// -- end extract COLUMNS of A(p,p) of block i --
	// ----------------------------------------------
	// ----------------------------------------------
	

	// ----------------------------------------------
	// ----------------------------------------------
        // --- extract columns of A(p,p)^T of block i ---
	// pointer to the current column inside diagonal block i
        pD=AvalD[i_bldl];
	// pointer to the current strict lower triangular block
	pL=AvalL;
	// reset j
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 start
        j=r; s=r+bi;
	// loop through the rows r:s-1 inside block column i and search
	// for the required columns
	for (k=0; k<bi; k++,pD+=bi,pL+=BL->nblockrow[i_bldl]) {
	    jj=p[j]; // shortcut
	    // check ROW j of A(p,p)	    
	    // since A is stored by columns we access row j via the linked list
	    jt=Ahead[j];
	    while (jt>=0) {
	          // printf("scanning column %ld\n",jt);
	          jjt=p[jt]; // shortcut
		  // check COLUMN jt of A(p,p), its leading entry must be j
	          m=Afirst[jt];
		  // row index tt of A(tt,p[jt])
		  tt=A->rowind[jjt][m];
		  // tt=q[t], we must have j=t, i.e. we are considerung A(p[j],p[jt])
		  t=invq[tt];
	          // printf("checking A(p(%ld),p(%ld))\n",t,jt);
		  
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
		     rval=SL[jjt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		     val*=rval;
#else
		     val.r*=rval;
		     val.i*=rval;
#endif
		  } // end if SL!=0
		  
		  // jt is inside the diagonal block
		  if (r<=jt && jt<s) {
		     // insert conjugate complex block-diagonal entry
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		     pD[jt-r]=val;
#else
		     pD[jt-r].r= val.r;
		     pD[jt-r].i=-val.i;
#endif
		  } // end if r<=jt<s
		  // jt is located below the diagonal block
		  else if (jt>=s) {
		     // check whether this is a new row or not
		     tt=idxposL[jt];
		     // new nonzero row index found in this block row
		     if (!tt) {
		        idxL[cntL]=jt;      // bookmark index and ...		
			idxposL[jt]=++cntL; // its physical position, shifted by 1
			tt=cntL;	    // return its assigned location       
			// if (nrowind>=nrAL) {
			//    nrAL+=10;
			//    Arowind=(integer *)realloc(Arowind,(size_t)nrAL*sizeof(integer));
			// } // end if
			// bookmark index for later purposes
		        // Arowind[nrowind++]=jt;
			nrowind++;
			// nrowind and cntL are identical; similarly, idxL and Arowind coincide
		     } // end if
		     // whoops, our estimate about the number of nonzeros in
		     // the current block column was too small
		     if (nrowind>BL->nblockrow[i]) {
		        // printf("whoops, step %ld, A^T, realloc\n",i);
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
			    for (mm=0; mm<(BL->nblockrow[i]-1); mm++,pU--) {
			        pL[jj]=*pL;
			    } // end for mm
			} // end for jj
			// set additional element in column 0 to 0.0
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			pL[1]=0.0;
#else
			pL[1].r=pL[1].i=0.0;
#endif
			// re-adjust pointer w.r.t. new number of rows
			pL=AvalL+k*BL->nblockrow[i];
		     } // end if
		     // insert conjugate complex sub-diagonal entry
		     tt-=bi; // remember that idxL, idxposL include diagonal block
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		     pL[tt-1]=val;
#else
		     pL[tt-1].r= val.r;
		     pL[tt-1].i=-val.i;
#endif
		  } // end else-if jt>=s      
	          jt=Alist[jt];
	    } // end while jt>=0
	    // ---------------------------------------------------------------
	    // ----- update linked list for the next row  greater than j -----
	    // beginning of the linked list
	    jt=Ahead[j];
	    // while list is non-empty
	    while (jt>=0) {
	          // position of the current leading entry in column jt referring to A(p[j:end],p[jt])
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
	      printf("%12ld",Ahead[m]);
	    printf("\n");
	    for (m=0; m<n; m++)
	      printf("%12ld",Alist[m]);
	    printf("\n");
	    for (m=0; m<n; m++)
	      printf("%12ld",Afirst[m]);
	    printf("\n");
	    */
	    
	    // advance row counter
	    j++;
	} // end for k

	// number of off-diagonal rows
	// nrowind=cntL-bi;
	/*
	if (nrowind!=cntL-bi) {
	   printf("step %ld, nrowind=%ld, cntL=%ld, bi=%ld\n",
		  i_bldl,nrowind,cntL,bi);
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

	// sort AvalL such that its row indices are sorted
	// in increasing order
	// QSORTGNL(AvalL,Arowind,idxbuff,&nrowind,&bi);
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
	
	// for testing, print block structure as extracted from A
#ifdef PRINT_INFO
	printf("block %3ld of A\n",i);fflush(stdout);
	printf("diagonal block of A\n            ");
	for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	    printf("%12ld",BiD->colind[i_bldl][m]);
	printf("\n");fflush(stdout);
	for (k=0; k<BiD->nblockcol[i_bldl]; k++) {
	    printf("%12ld",BiD->colind[i_bldl][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	        printf("%12.4le",AvalD[i_bldl][m*BiD->nblockcol[i_bldl]+k]);
	    printf("\n");fflush(stdout);
#else	    
	    for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	        printf("%12.4le",AvalD[i_bldl][m*BiD->nblockcol[i_bldl]+k].r);
	    printf("\n            ");
	    for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	        printf("%12.4le",AvalD[i_bldl][m*BiD->nblockcol[i_bldl]+k].i);
	    printf("\n");fflush(stdout);
#endif
	}
	printf("\n");

	printf("sub-diagonal block of A\n            ");
	if (nrowind>0) {
	   for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	       printf("%12ld",BiD->colind[i_bldl][m]);
	   printf("\n");
	   for (k=0; k<nrowind; k++) {
	       printf("%12ld",idxL[bi+k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
	           printf("%12.4le",AvalL[m*nrowind+k]);
#else
	       for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
		   printf("%12.4le",AvalL[m*nrowind+k].r);
	       printf("\n            ");
	       for (m=0; m<BiD->nblockcol[i_bldl]; m++) 
		   printf("%12.4le",AvalL[m*nrowind+k].i);
#endif
	       printf("\n");
	   } // end for k
	}
	printf("\n");
#endif
        // - end extract COLUMNS of A(p,p)^T of block i -
	// ----------------------------------------------
	// ----------------------------------------------


#ifdef _PROFILING_
	time_extract_block+=omp_get_wtime()-timeBegin;
#endif
	



	// now that block column i has been extracted into two blocks
	// we can compute block column i_bldl of L
	// To do so, we use linked list for rows of BL associated with the
	// block columns of BL 
	

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// ----------------------------------------------------------------------
	// compute block column i_bldl of BL using linked list for the rows of BL
	// r,s refer to the column/row indices where diagonal block
	// i and i+1 of A start
        r=BiD->colind[i_bldl][0]; s=r+bi;

	
	// first pass: set up nonzero row structure of D+L before dropping
	// store indices of the diagonal block of A
	// for (m=r; m<s; m++) {
	//    idxL[cntL]=m;
	//    idxposL[m]=++cntL;
	// } // end for m
	// store indices of the sub-diagonal block of A
	// for (m=0; m<nrowind; m++) {
	//    t=Arowind[m];
	//    idxL[cntL]=t;
	//    idxposL[t]=++cntL;
	// } // end for m

	
	// beginning of the linked list for block column i
	// this information is obtained from the sparsity pattern of the rows of L
	it=BLhead[i_bldl];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning block column %ld of BL\n",it); fflush(stdout);
#endif
	      // check block COLUMN it of BL
	      rowind=BL->rowind[it];
	      nrow=BL->nblockrow[it];
	      if (cntL>bi) {
		 // idxL includes more indices than just the diagonal block
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
#ifdef PRINT_INFO
			printf("new row index %ld added\n",t); 
#endif
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
			  idxposL[t]=i2--;
			  idxL[i2]=t;
			  ii--;
			  t=idxU[ii]; 
		       }
		       else { // tt>t
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

	      
	      /*
	      if (0<cntL_old && cntL_old<cntL) {
		 ii=0; jj=cntL_old;
		 while (ii<cntL_old && jj<cntL) {
		       t=idxL[ii]; tt=idxL[jj];
		       if (t<tt) {
			  idxU[cntU++]=t;
			  ii++;
		       }
		       else { // tt<t, the case tt=t is avoided by flagging idxposL
			  idxU[cntU++]=tt;
			  jj++;
		       }
		 } // end while
		 while (ii<cntL_old) {
		       t=idxL[ii];
		       idxU[cntU++]=t;
		       ii++;
		 } // end while
		 while (jj<cntL) {
		       tt=idxL[jj];
		       idxU[cntU++]=tt;
		       jj++;
		 } // end while
		 cntU=0;
		 // write the sorted list back to idxL
		 for (m=0; m<cntL; ) {
		     t=idxU[m];
		     idxL[m]=t;
		     idxposL[t]=++m;
		 } // end for m
	      } // end if 0<cntL_old<cntL
	      */
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		if (idxposU[m]) {
		  printf("5.idxposU[%ld]=%ld\n",m,idxposU[m]);
		  fflush(stdout);
		} // end if
	      } // end for m	
#endif
	      // advance to next required block column of BL as requested by BL^T
	      it=BLlist[it];
	} // end while
	// end first pass
#ifdef PRINT_INFO
	printf("total nonzero row indices column %ld of BL\n",i); fflush(stdout);
	for (m=0; m<cntL; m++) {
  	    printf("%12ld",idxL[m]);
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
	pb=buff; pD=AvalD[i_bldl];
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
	time_LDL_update_pass_1+=omp_get_wtime()-timeBegin;
#endif
	
	// second pass: compute block column i of BL using level-3-BLAS
	// beginning of the linked list for the rows of BL which 
	// reflects the requested block columns of BL
	it=BLhead[i_bldl];
	// while list is non-empty
	while (it>=0) {
#ifdef PRINT_INFO
	      printf("scanning again block column %ld of BL\n",it); fflush(stdout);
#endif
	      // level-3-BLAS update "C+=AB^T"

#ifdef _PROFILING_
	      timeBegin = omp_get_wtime();
#endif
	      // step 1: rows of block column "it" of BL refer to rows of "B"
	      //         These can be left in-place for the level-3-BLAS update
	      // read the associated rows from BL for "B"
	      jj=BLfirst[it];        // start of row r and higher
	      ldb=BL->nblockrow[it]; // leading dimension (total number of rows)
	      flag=-1;
	      rowind=BL->rowind[it]; // shortcut
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
#ifdef PRINT_INFO
	      printf("column indices required to update with column %ld of BL\n",it); fflush(stdout);
	      for (m=0; m<cntU; m++) {
		  printf("%12ld",idxU[m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef PRINT_INFO
	      printf("row indices required to update with column %ld of BL\n",it); fflush(stdout);
	      for (m=BLfirst[it]; m<BL->nblockrow[it]; m++) {
		  printf("%12ld",BL->rowind[it][m]);
	      } // end for m
	      printf("\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	      time_LDL_update_pass_2+=omp_get_wtime()-timeBegin;
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
		 time_LDL_update_pass_2+=omp_get_wtime()-timeBegin;
#endif
		 if (cntU==1 && l==1) {
#ifdef PRINT_INFO
		    printf("scalar update with block column %ld\n",it);fflush(stdout);
#endif
#ifdef _PROFILING_
		    timeBegin=omp_get_wtime();
#endif
		    
		    // extract associated diagonal entry
		    valD=AvalD[it][0];
		    // extract column index t
		    t=idxU[0];
		    // scatter to column t of buff
		    pb=buff+cntL*(t-r)-1;  // associated column of buff,shifted
		    pL=BL->valE[it] +ii; // numerical values of column BL{it}
		    // pU=BL->valE[it] +BLfirst[it]; // numerical value BL{it}^T
		    pU=pL; // numerical value BL{it}^T
		    
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
		    // scale BL{it}^T by the negative diagonal entry
		    val=-valD*(*pU);
#else
		    // scale CONJG(BL{it}^T) by the negative diagonal entry
		    val.r=-valD.r*pU->r-valD.i*pU->i;
		    val.i= valD.r*pU->i-valD.i*pU->r;
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
			// downdate buff - BL{it}*CONJG(BL{it}^T)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			pb[i2]+=(*pL++)*val;
#else
			// now update
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
		    // pU=BL->valE[it] +BLfirst[it]; // numerical value BL{it}^T
		    pU=BL->valE[it] +ii; // numerical value BL{it}^T
		    if (l==1) {
#ifdef PRINT_INFO
		       printf("small block rank-1 update with block column %ld\n",it);fflush(stdout);
#endif
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
		       // extract associated diagonal entry
		       valD=AvalD[it][0];
 		       for (m=0; m<cntU; m++,pU++) {
			   // numerical values of column BL{it}
			   pL=BL->valE[it]  +ii; 
			   // extract column index t
			   t=idxU[m];
			   // directly scatter to column t of buff
			   pb=buff+cntL*(t-r)-1; // associated column of buff, shifted
			   
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   // scale BL{it}^T by the negative diagonal entry
			   val=-valD*(*pU);
#else
			   // scale CONJG(BL{it}^T) by the negative diagonal entry
			   val.r=-valD.r*pU->r-valD.i*pU->i;
			   val.i= valD.r*pU->i-valD.i*pU->r;
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
			       // downdate buff - BL{it}*BL{it}^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
			       pb[i2]+=(*pL++)*val;
#else
			       // now update
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
#ifdef PRINT_INFO
		       printf("small block level-1 and level-2 update with block column %ld\n",it);fflush(stdout);
#endif
#ifdef _PROFILING_
		       timeBegin=omp_get_wtime();
#endif
		       for (m=0; m<cntU; m++,pU++) {
			   // numerical values of column BL{it}
			   pL=BL->valE[it]  +ii; 
			   // extract column index t
			   t=idxU[m];
			   // directly scatter to column t of buff
			   pb=buff+cntL*(t-r)-1; // associated column of buff, shifted

			   // scale by the diagonal block
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   for (mm=0; mm<l; mm++) 
			       buffD[mm]=0.0;
			   alpha=1.0; beta=0.0;
#else
			   for (mm=0; mm<l; mm++) 
			       buffD[mm].r=buffD[mm].i=0.0;
			   alpha.r=1.0; beta.r=0.0; alpha.i=beta.i=0.0;
#endif
#ifdef PRINT_INFO1
			   printf("diagonal block D(%ld,%ld)\n",it,it);
			   for (i2=0; i2<l; i2++) {
			       for (mm=0; mm<l; mm++)
				   printf("%12.4le",AvalD[it][l*mm+i2]);
			       printf("\n");fflush(stdout);
			   }
			   printf("\n");fflush(stdout);
			   printf("BL{%ld}(%ld,:)^T\n",it,idxU[m]);fflush(stdout);
			   for (mm=0; mm<l; mm++)
			       printf("%12.4le\n",pU[mm*ldb]);
			   printf("\n");fflush(stdout);
				   
#endif
			   mm=1;
			   // buffD = 1*D_{it,it}*BL{it}(m,:)^T+0*buffD
			   transa="t";
			   GEMV(transa,&l,&l,&alpha, AvalD[it],&l, pU,&ldb,
				&beta, buffD,&mm, 1);
#ifdef PRINT_INFO1
			   printf("diagonal block D(%ld,%ld)*BL{%ld}(%ld,:)^T\n",it,it,it,idxU[m]);fflush(stdout);
			   for (mm=0; mm<l; mm++)
			       printf("%12.4le\n",buffD[mm]);
			   printf("\n");fflush(stdout);
			   mm=1;				   
#endif
			   for (k=0; k<ldc; k++) {
			       // position i2 in buff
			       i2=idxbuff[k];
			       // downdate buff - BL{it}*(D{it,it}*BL{it}^T)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			       pb[i2]-=DISTDOT(&l, buffD,&mm, pL,&lda);
#else
#ifdef _USE_MKL_
			       DISTDOT(&val, &l, buffD,&mm, pL,&lda);
#else
			       val=DISTDOT(&l, buffD,&mm, pL,&lda);
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
#ifdef PRINT_INFO
		    printf("level-3 BLAS update with block column %ld\n",it);fflush(stdout);
#endif
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
			       // position i2 in buff
			       i2=idxbuff[k];
			       *pC++=pb[i2];
			   } // end for k
#endif //-else _USE_MKL_
			
		       } // end for m
		    } // end if NOT skip_gthr_sctr
		    
		    // now perform level-3-BLAS updates
		    
		    // 1. buffD = 1*D*CONJG(B^T)+0*buffD
		    if (l*cntU>nbuffD) {
		       nbuffD=MAX(nbuffD+100,l*cntU);
		       buffD=(FLOAT *) realloc(buffD,nbuffD*sizeof(FLOAT));
		    } // end if
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    for (mm=0; mm<l*cntU; mm++) 	        
		        buffD[mm]=0.0;
		    alpha=1.0;
		    beta =0.0;
		    
#else
		    for (mm=0; mm<l*cntU; mm++) 	        
		        buffD[mm].r=buffD[mm].i=0.0;
		    alpha.r= 1.0; alpha.i=0.0;
		    beta.r = 0.0; beta.i =0.0;
#endif
		    transa="n"; transb="c";
		    GEMM(transa,transb, &l,&cntU,&l, &alpha,
			 AvalD[it],&l,
			 BL->valE[it]+BLfirst[it],&ldb, &beta, buffD,&l, 1,1);
		    
		    // 2. C = -1*A*buffD+1*C
		    transa="n"; transb="n";
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
			    buffD,&l, &beta,
			    buff + idxbuff[0]-1 + cntL*(idxU[0]-r),&cntL, 1,1);
		    else // level-3 BLAS uses a cached copy
		       GEMM(transa,transb, &ldc,&cntU,&l, &alpha,
			    BL->valE[it]+BLfirst[it],&lda,
			    buffD,&l, &beta,
			    buffC,&ldc, 1,1);
	      
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
			       // position of i2 in buff
			       i2=idxbuff[k];
			       pb[i2]=*pC++;
			   } // end for k
#endif //-else _USE_MKL_
			
		       } // end for m
		    } // end if NOT skip_gthr_sctr
#ifdef _PROFILING_
		    time_level_3_blas_update+=omp_get_wtime()-timeBegin;
#endif
		 } // end if-elseif-else cntU=1 and l=1
	      } // end if ldc> and cntU>0 and l>0
		 
#ifdef _PROFILING_
	      timeBegin=omp_get_wtime();
#endif
	      // update BLfirst, this automatically excludes the diagonal
	      // block when block ROW i is updated
	      BLfirst[it]=jj;

	      // reset idxposU
	      for (m=0; m<cntU; m++) {
		  t=idxU[m];
		  idxposU[t]=0;
	      } // end for m
	      cntU=0;
#ifdef PRINT_CHECK
	      for (m=0; m<n; m++) {
		if (idxposU[m]) {
		  printf("6.idxposU[%ld]=%ld\n",m,idxposU[m]);
		  fflush(stdout);
		} // end if
	      } // end for m	
#endif
	      
	      // advance to the next block column
	      it=BLlist[it];
#ifdef _PROFILING_
	      time_LDL_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	} // end while
	// end second pass
#ifdef PRINT_INFO
	printf("computed block column %ld of BL+D\n",i_bldl);fflush(stdout);
	printf("        ");
	for (m=r; m<s; m++)
	    printf("%12ld",m);
	printf("\n");fflush(stdout);
	for (ii=0; ii<cntL; ii++) {
	    printf("%12ld",idxL[ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",buff[m*cntL+ii]);
	    }
	    printf("\n");fflush(stdout);
#else
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",buff[m*cntL+ii].r);
	    }
	    printf("\n            ");fflush(stdout);
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",buff[m*cntL+ii].i);
	    }
	    printf("\n");fflush(stdout);
#endif
	}
	printf("\n");fflush(stdout);
#endif

#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif

	// force diagonal entries to be real-valued
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	pb=buff;
	for (m=0; m<bi; m++,pb+=cntL+1)
	    pb->i=0.0;
#endif
	// for technical reasons, make a copy of the diagonal block
	// (stored in the begining of buff) to AvalD[i_bldl] before inversion
	
	pb=buff; pD=AvalD[i_bldl];
	for (m=0; m<bi; m++,pb+=cntL,pD+=bi) {
	    // re-insert diagonal entries
	    if (perturbation)
	       pb[m]=buffDgl[idxL[m]];
	    // copy columns r,...,s-1 of the diagonal block
	    memcpy(pD,pb,(size_t)bi*sizeof(FLOAT));
	} // end for m
#ifdef PRINT_INFO
	printf("copy of the original diagonal block BD{%ld}\n",i_bldl);fflush(stdout);
	printf("        ");
	for (m=r; m<s; m++)
	    printf("%12ld",m);
	printf("\n");fflush(stdout);
	for (ii=0; ii<bi; ii++) {
	    printf("%12ld",idxL[ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",AvalD[i_bldl][m*bi+ii]);
	    }
	    printf("\n");fflush(stdout);
#else
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",AvalD[i_bldl][m*bi+ii].r);
	    }
	    printf("\n            ");
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",AvalD[i_bldl][m*bi+ii].i);
	    }
	    printf("\n");fflush(stdout);
#endif
	}
	printf("\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	
	// scalar column, skip LAPACK
	if (bi==1) {
	  
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // the diagonal entry is buff[0]
  	   rval=FABS2(buff[0]); // |diagonal_entry|^2
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	   if (buff[0]<=0.0)
#else
	   if (buff[0].r<=0.0)
#endif
	   {
	      ierr=1;
	      printf("BICHOL encountered zero or negative pivot\n");

	      // sort original matrix back such that row indices and their
	      // values are taken in increasing order
	      for (j=0; j<n; j++)
		  QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
	      
	      // printf("0.0\n"); fflush(stdout);
	      idxL--;idxL=FREE(idxL);
	      // printf("0.01\n"); fflush(stdout);
	      idxU--;idxU=FREE(idxU);
	      // printf("0.02\n"); fflush(stdout);
	      idxbuff =FREE(idxbuff);
	      // printf("0.03\n"); fflush(stdout);
	      idxposL =FREE(idxposL);
	      // printf("0.04\n"); fflush(stdout);
	      idxposU =FREE(idxposU);
	      // printf("0.05\n"); fflush(stdout);
	      invblock=FREE(invblock);
	      // printf("0.06\n"); fflush(stdout);
	      Ahead   =FREE(Ahead);
	      // printf("0.07\n"); fflush(stdout);
	      Alist   =FREE(Alist);
	      // printf("0.08\n"); fflush(stdout);
	      Afirst  =FREE(Afirst);
	      // printf("0.09\n"); fflush(stdout);
	      BLhead  =FREE(BLhead);
	      // printf("0.010\n"); fflush(stdout);
	      BLlist  =FREE(BLlist);
	      // printf("0.011n"); fflush(stdout);
	      BLfirst =FREE(BLfirst);
	      // printf("0.012\n"); fflush(stdout);
	      AvalL   =FREE(AvalL);
	      // printf("0.013\n"); fflush(stdout);
	      Arowind =FREE(Arowind);
	      // printf("0.014\n"); fflush(stdout);
	      buffC   =FREE(buffC);
	      // printf("0.015\n"); fflush(stdout);
	      buff    =FREE(buff);
	      // printf("0.016\n"); fflush(stdout);
	      buffD   =FREE(buffD);
	      // printf("0.1\n"); fflush(stdout);

	      for (i=0; i<nblocks; i++) {
		  BL->colind[i] =FREE(BL->colind[i] );
		  BL->rowind[i] =FREE(BL->rowind[i] );
		  BiD->colind[i]=FREE(BiD->colind[i]);
		  BL->valE[i]   =FREE(BL->valE[i]   );
		  BiD->valD[i]  =FREE(BiD->valD[i]  );
		  AvalD[i]      =FREE(AvalD[i]      );
	      }
	      AvalD   =FREE(AvalD);
	      // printf("0.2\n"); fflush(stdout);
	      
	      BL->nblockcol =FREE(BL->nblockcol);
	      BL->nblockrow =FREE(BL->nblockrow);
	      BL->colind    =FREE(BL->colind);
	      BL->rowind    =FREE(BL->rowind);
	      BL->valE      =FREE(BL->valE);

	      BiD->nblockcol=FREE(BiD->nblockcol);
	      BiD->colind   =FREE(BiD->colind   );
	      BiD->valD     =FREE(BiD->valD     );
	      
	      BL->nc=BL->nr=BiD->nc=BiD->nr=0;
	      BL->nblocks=BiD->nblocks=0;
	      // printf("0.3\n"); fflush(stdout);
	      
	      if (flag_ildl1t) {
		 free(blocksize);
		 blocksize=myblocksize;
		 nblocks  =mynblocks;
	      } // end if flag_ildl1t
	      // printf("0.4\n"); fflush(stdout);
	      
	      if (perturbation)
		 buffDgl=FREE(buffDgl);
	      // printf("0.5\n"); fflush(stdout);
	      
	      return (ierr);
	   }
	   

	   
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
	   BiD->valD[i_bldl]=(FLOAT *)malloc(sizeof(FLOAT));
	   BiD->valD[i_bldl][0]=buff[0];
#ifdef _PROFILING_
	   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // apply dropping to the strict lower triangular block of BL
	   ldc=cntL-1; // number of rows before dropping
	   srval=sqrt(rval);
	   rval=srval*droptol; // |diagonal_entry|*droptol
	   pb=buff+1; // start strict lower triangular part
	   l=0;
	   if (perturbation) {
	      // some inits for Ajiz&Jennings
	      srval=sqrt(srval);  // sqrt(diagonal_entry)
	      e=0.0; // additional lumping for the diagonal entry
	      ierr=0;
	      for (ii=0; ii<ldc; ii++,pb++) {
		  ab=FABS(*pb);
		  if (ab>=rval)
		     l++;
		  else {
		     // add weight to the remaining diagonal entries
		     jj=idxL[ii+1]; // associated index shifted by 1 to skip the diagonal index
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     sb=sqrt(buffDgl[jj]);
		     if (perturbation>0)
		        buffDgl[jj]+=(sb/srval)*ab;
		     else
		        buffDgl[jj]-=(sb/srval)*ab;
		     if (buffDgl[jj]<=0.0)
		        ierr=1;
#else
		     sb=sqrt(buffDgl[jj].r);
		     if (perturbation>0)
		        buffDgl[jj].r+=(sb/srval)*ab;
		     else
		        buffDgl[jj].r-=(sb/srval)*ab;
		     if (buffDgl[jj].r<=0.0)
		        ierr=1;
#endif
		     // accumulate weight for the current diagonal entry
		     if (perturbation>0)
		        e=e+ab/sb;
		     else
		        e=e-ab/sb;
		  } // end if-else
	      } // end for ii
	      
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      // update the copy of the diagonal entry
	      AvalD[i_bldl][0]+=srval*e;
	      // recompute inverse
	      buff[0]=1.0/AvalD[i_bldl][0];
	      if (AvalD[i_bldl][0]<=0.0)
		 ierr=1;
#else
	      // update the copy of the diagonal entry
	      AvalD[i_bldl][0].r+=srval*e;
	      // recompute inverse
	      buff[0].r=1.0/AvalD[i_bldl][0].r;
	      if (AvalD[i_bldl][0].r<=0.0)
		 ierr=1;
#endif
	      if (ierr) {
		 ierr=1;
		 printf("BICHOL encountered zero or negative diagonal part\n");


		 // sort original matrix back such that row indices and their
		 // values are taken in increasing order
		 for (j=0; j<n; j++)
		     QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
		 
		 // printf("2.0\n"); fflush(stdout);
		 idxL--;idxL=FREE(idxL);
		 // printf("2.01\n"); fflush(stdout);
		 idxU--;idxU=FREE(idxU);
		 // printf("2.02\n"); fflush(stdout);
		 idxbuff =FREE(idxbuff);
		 // printf("2.03\n"); fflush(stdout);
		 idxposL =FREE(idxposL);
		 // printf("2.04\n"); fflush(stdout);
		 idxposU =FREE(idxposU);
		 // printf("2.05\n"); fflush(stdout);
		 invblock=FREE(invblock);
		 // printf("2.06\n"); fflush(stdout);
		 Ahead   =FREE(Ahead);
		 // printf("2.07\n"); fflush(stdout);
		 Alist   =FREE(Alist);
		 // printf("2.08\n"); fflush(stdout);
		 Afirst  =FREE(Afirst);
		 // printf("2.09\n"); fflush(stdout);
		 BLhead  =FREE(BLhead);
		 // printf("2.010\n"); fflush(stdout);
		 BLlist  =FREE(BLlist);
		 // printf("2.011\n"); fflush(stdout);
		 BLfirst =FREE(BLfirst);
		 // printf("2.012\n"); fflush(stdout);
		 AvalL   =FREE(AvalL);
		 // printf("2.013\n"); fflush(stdout);
		 Arowind =FREE(Arowind);
		 // printf("2.014\n"); fflush(stdout);
		 buffC   =FREE(buffC);
		 // printf("2.015\n"); fflush(stdout);
		 buff    =FREE(buff);
		 // printf("2.016\n"); fflush(stdout);
		 buffD   =FREE(buffD);
		 // printf("2.1\n"); fflush(stdout);

		 for (ii=0; ii<nblocks; ii++) {
		     BL->colind[ii] =FREE(BL->colind[ii] );
		     BL->rowind[ii] =FREE(BL->rowind[ii] );
		     BiD->colind[ii]=FREE(BiD->colind[ii]);
		     BL->valE[ii]   =FREE(BL->valE[ii]   );
		     BiD->valD[ii]  =FREE(BiD->valD[ii]  );
		     AvalD[ii]      =FREE(AvalD[ii]      );
		 }
		 AvalD   =FREE(AvalD);
		 // printf("2.2\n"); fflush(stdout);
	      
		 BL->nblockcol =FREE(BL->nblockcol);
		 BL->nblockrow =FREE(BL->nblockrow);
		 BL->colind    =FREE(BL->colind);
		 BL->rowind    =FREE(BL->rowind);
		 BL->valE      =FREE(BL->valE);
		 
		 BiD->nblockcol=FREE(BiD->nblockcol);
		 BiD->colind   =FREE(BiD->colind   );
		 BiD->valD     =FREE(BiD->valD     );
		 
		 BL->nc=BL->nr=BiD->nc=BiD->nr=0;
		 BL->nblocks=BiD->nblocks=0;
		 // printf("2.3\n"); fflush(stdout);
		 
		 if (flag_ildl1t) {
		    free(blocksize);
		    blocksize=myblocksize;
		    nblocks  =mynblocks;
		 } // end if flag_ildl1t
		 // printf("2.4\n"); fflush(stdout);
	      
		 if (perturbation)
		    buffDgl=FREE(buffDgl);
		 // printf("2.5\n"); fflush(stdout);
	      
		 return (ierr);
	      } // end if ierr
	   } // end if perturbation
	   else { // no lumping
	      // skip dropping check if useless
	      if (droptol<=0.0)
		 l=ldc;
	      else
		 for (ii=0; ii<ldc; ii++,pb++)
		     if (FABS(*pb)>=rval)
		        l++;
	   } // end if-else perturbation
	   
	   // now we know the fill
	   BL->nblockrow[i_bldl]=l;
	   BL->valE[i_bldl]=(FLOAT *)malloc((size_t)l*sizeof(FLOAT));
	   BL->rowind[i_bldl]=(integer *)malloc((size_t)l*sizeof(integer));
	   	   
	   // multiply strict lower triangular part by the inverse diagonal entry
	   // and copy it to BL
	   pb=buff+1;
	   pL=BL->valE[i_bldl];
	   pi=BL->rowind[i_bldl];
	   if (perturbation) {
	      for (ii=0; ii<ldc; ii++,pb++)
	           if (FABS(*pb)>=rval) {
		      jj=idxL[ii+1]; // shift idxL by 1 because of the diagonal index
		      *pi++=jj;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		  
		      d=(*pb)*(*buff);
		      // downdate associated diagonal entry as usual
		      buffDgl[jj]-=d*(*pb);
		      // copy buff[ii]/D to L
		      *pL++=d;
#else
		      // downdate associated diagonal entry as usual
		      d=pb->r*buff->r;
		      buffDgl[jj].r-=d*pb->r;
		      d=pb->i*buff->r;
		      buffDgl[jj].r-=d*pb->i;
		      // copy buff[ii]/D to L
		      pL->r=pb->r*buff->r-pb->i*buff->i;
		      pL->i=pb->r*buff->i+pb->i*buff->r;
		      pL++;
#endif		  
		   } // end if |*pb|>=rval
#ifdef _PROFILING_
	      time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	   } // end if perturbation
	   else { // no lumping
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
		  } // end if |*pb|>=rval
#ifdef _PROFILING_
	      time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	   } // end if-else perturbation
        } // end if bi=1
	else { // bi>1 block case
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // factorize and invert diagonal block using LAPACK
	   POTRF("l",&bi, buff,&cntL, &ierr,1);
	   if (ierr) {
	      if (ierr<0)
		 printf("LAPACK's POTRF: %ld-th argument had an illegal value\n", -ierr);
	      else
		 printf("LAPACK's POTRF routine encountered zero column in step %ld\n", ierr);

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
	      AvalL   =FREE(AvalL);
	      Arowind =FREE(Arowind);
	      buffC   =FREE(buffC);
	      buff    =FREE(buff);
	      buffD   =FREE(buffD);

	      for (i=0; i<nblocks; i++) {
		  BL->colind[i] =FREE(BL->colind[i] );
		  BL->rowind[i] =FREE(BL->rowind[i] );
		  BiD->colind[i]=FREE(BiD->colind[i]);
		  BL->valE[i]   =FREE(BL->valE[i]   );
		  BiD->valD[i]  =FREE(BiD->valD[i]  );
		  AvalD[i]      =FREE(AvalD[i]      );
	      }
	      AvalD   =FREE(AvalD);
	      
	      BL->nblockcol =FREE(BL->nblockcol);
	      BL->nblockrow =FREE(BL->nblockrow);
	      BL->colind    =FREE(BL->colind);
	      BL->rowind    =FREE(BL->rowind);
	      BL->valE      =FREE(BL->valE);

	      BiD->nblockcol=FREE(BiD->nblockcol);
	      BiD->colind   =FREE(BiD->colind   );
	      BiD->valD     =FREE(BiD->valD     );
	      
	      BL->nc=BL->nr=BiD->nc=BiD->nr=0;
	      BL->nblocks=BiD->nblocks=0;
	      
	      if (flag_ildl1t) {
		 free(blocksize);
		 blocksize=myblocksize;
		 nblocks  =mynblocks;
	      } // end if flag_ildl1t
	      
	      if (perturbation)
		 buffDgl=FREE(buffDgl);
	      
	      return (ierr);
	   } // end if ierr


	   // compute local determinant
	   locdet.r=0.0;
	   locdet.i=0.0;
	   ii=0;
	   for (jj=0; jj<bi; jj++,ii+=cntL) {
	       // 1x1 pivot
	       ajj=buff[ii+jj];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       // mexPrintf("ajj=%12.4le\n",ajj);
	       if (ajj>=0.0)
		  locdet.r+=2.0*log(ajj);
	       else { // log of negative numbers is complex (main branch)
		  locdet.r+=2.0*log(-ajj); 
		  locdet.i+=M_PI;
	       }
#else
	       // mexPrintf("ajj=%12.4le+i%12.4le\n",ajj.r,ajj.i);
	       // complex logarithm (main branch)
	       rval=FABS(ajj);
	       // 1. real part
	       locdet.r+=2.0*log(rval);
	       // 1. imaginary part
	       // compute argument of ajj
	       if (ajj.i>=0) 
		  rval= 2.0*acos(ajj.r/rval);
	       else
		  rval=-2.0*acos(ajj.r/rval);
	       locdet.i+=rval;	       
#endif
	   } // end for jj

	   
	   if (invert_blocks) {
	      POTRI("l",&bi, buff,&cntL, &ierr,1);
	      if (ierr) {
		 if (ierr<0)
		    printf("LAPACK's POTRI: %ld-th argument had an illegal value\n", -ierr);
		 else
		    printf("LAPACK's POTRI routine encountered zero column in step %ld\n", ierr);

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
		 AvalL   =FREE(AvalL);
		 Arowind =FREE(Arowind);
		 buffC   =FREE(buffC);
		 buff    =FREE(buff);
		 buffD   =FREE(buffD);

		 for (i=0; i<nblocks; i++) {
		     BL->colind[i] =FREE(BL->colind[i] );
		     BL->rowind[i] =FREE(BL->rowind[i] );
		     BiD->colind[i]=FREE(BiD->colind[i]);
		     BL->valE[i]   =FREE(BL->valE[i]   );
		     BiD->valD[i]  =FREE(BiD->valD[i]  );
		     AvalD[i]      =FREE(AvalD[i]      );
		 }
		 AvalD   =FREE(AvalD);
		 
		 BL->nblockcol =FREE(BL->nblockcol);
		 BL->nblockrow =FREE(BL->nblockrow);
		 BL->colind    =FREE(BL->colind);
		 BL->rowind    =FREE(BL->rowind);
		 BL->valE      =FREE(BL->valE);
		 
		 BiD->nblockcol=FREE(BiD->nblockcol);
		 BiD->colind   =FREE(BiD->colind   );
		 BiD->valD     =FREE(BiD->valD     );
	      
		 BL->nc=BL->nr=BiD->nc=BiD->nr=0;
		 BL->nblocks=BiD->nblocks=0;
	      
		 if (flag_ildl1t) {
		    free(blocksize);
		    blocksize=myblocksize;
		    nblocks  =mynblocks;
		 } // end if flag_ildl1t
	      
		 if (perturbation)
		    buffDgl=FREE(buffDgl);
	      
		 return (ierr);
	      } // end if ierr

	      // copy strict lower triangular part to strict upper triangular part
	      pb=buff+1;
	      mm=1;
	      for (ii=1; ii<bi; pb++,ii++) {
		  pD=buff+cntL*ii;
		  COPY(&ii, pb,&cntL, pD,&mm);
		  // conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
		  for (m=0; m<ii; pD++,m++)
		      pD->i=-pD->i;
#endif
	      } // end for ii
	   } // end if invert_blocks
	   else {
	      // nothing, no inversion, no copying
	   }
	   
	
	   
	   // extract (inverse) diagonal block to BiD
	   BiD->valD[i_bldl]=(FLOAT *)malloc((size_t)bi*bi*sizeof(FLOAT));
	   pD=BiD->valD[i_bldl];
	   pb=buff;
	   for (ii=0; ii<bi; ii++,pD+=bi,pb+=cntL)
	       memcpy(pD,pb,(size_t)bi*sizeof(FLOAT));
	   
#ifdef PRINT_INFO0
	   if (invert_blocks)
	      printf("computed inverse diagonal block %ld of BiD\n",i_bldl);
	   else
	      printf("factorized diagonal block %ld of BiD\n",i_bldl);
	   printf("        ");
	   for (m=0; m<bi; m++)
	       printf("%12ld",BiD->colind[i_bldl][m]);
	   printf("\n");fflush(stdout);
	   for (ii=0; ii<bi; ii++) {
	       printf("%12ld",BiD->colind[i_bldl][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	       for (m=0; m<bi; m++) {
		   printf("%12.4le",BiD->valD[i_bldl][m*bi+ii]);
	       }
	       printf("\n");fflush(stdout);
#else
	       for (m=0; m<bi; m++) {
	           printf("%12.4le",BiD->valD[i_bldl][m*bi+ii].r);
	       }
	       printf("\n            ");
	       for (m=0; m<bi; m++) {
	           printf("%12.4le",BiD->valD[i_bldl][m*bi+ii].i);
	       }
	       printf("\n");fflush(stdout);
#endif
	   }
	   printf("\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   // allocate memory for BL
	   ldc=cntL-bi; // number of rows before dropping
	   k=ldc*bi;    // total size
	   BL->valE[i_bldl]=(FLOAT *)malloc((size_t)k*sizeof(FLOAT));
	   if (invert_blocks) {
	      // multiply strict lower triangular block by the inverse diagonal block
	      // and copy it to BL, use level-3-BLAS 
	      // init with 0
	      pL=BL->valE[i_bldl];
	      for (ii=0; ii<k; ii++,pL++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  *pL=0.0;
#else
		  pL->r=pL->i=0.0;
#endif
	      } // end for ii
	      // C=BL(:,i), A=subdgl(buff)=subdgl([BiD(i,i);BL(:,i)]), B=BiD(i,i)
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
		      BiD->valD[i_bldl],&bi, &beta, BL->valE[i_bldl],&ldc,1,1);
	   } // end if invert_blocks
	   else { // Cholesky case
	      // forward/backward solve using the Cholesky factorization of the diagonal block
	      // copy subdiagonal part of buff to BL since solve works in-place with BL^H
	      m=bi*ldc;
	      if (m>nbuffC) {
		 nbuffC=ELBOW_BUFF*m+BUFF_EXT;
		 buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	      } // end if
	      // buffC = subdgl(buff)^H = subdgl([BiD(i,i);BL(:,i)])^H
	      pb=buff+bi;
	      pL=buffC;
	      i2=1;
	      for (ii=0; ii<ldc; ii++,pb++,pL+=bi) 
		  COPY(&bi,pb,&cntL,pL,&i2);
	      // conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	      pL=buffC;
	      for (ii=0; ii<m; ii++,pL++) 
		  pL->i=-pL->i;
#endif
	      // forward/backward solves BiD(i,i)*BiD(i,i)^H * BL^H  <- BL^H
	      //                    <=>  BiD(i,i)*BiD(i,i)^H * buffC <- buffC
	      POTRS("l", &bi, &ldc, BiD->valD[i_bldl],&bi, buffC,&bi, &ierr,1);
	      if (ierr<0) {
		 printf("LAPACK's POTRS: %ld-th argument had an illegal value\n", -ierr);
		 return (ierr);
	      }
	      // BL(:,i)=buffC^H
	      pb=buffC;
	      pL=BL->valE[i_bldl];
	      i2=1;
	      for (ii=0; ii<ldc; ii++,pb+=bi,pL++) 
		  COPY(&bi,pb,&i2,pL,&ldc);
	      // conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	      pL=BL->valE[i_bldl];
	      for (ii=0; ii<m; ii++,pL++) 
		  pL->i=-pL->i;
#endif
	   } // end if-else invert_blocks

	   
	   if (perturbation) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      transa="t"; transb="n";
	      alpha=1.0;
	      beta =1.0;
#else
	      transa="c"; transb="n";
	      alpha.r=1.0; alpha.i=0.0;
	      beta.r =1.0; beta.i =0.0;
#endif
	      // apply dropping to the strict lower triangular block of BL
	      pL=BL->valE[i_bldl];
	      pb=buff+bi; // original entries are still in "buff", skip dgl blk
	      cntU=0;
	      if (BL->nblockrow[i_bldl]<ldc) {
		 // temporarily advance idxL
		 idxL+=bi;
		 for (ii=0; ii<ldc; ii++,pL++,pb++) {
		     jj=I_AMAX(&bi,pL,&ldc)-1;
		     k=jj*ldc;
		     // keep large entries
		     if (FABS(pL[k])>=droptol) {
		     } // if
		     else { // drop and compensate small entries
		        // get true row index
		        jj=idxL[ii];
			// and keep it in mind
			idxU[cntU++]=jj;
			// keep in mind which index position was chosen, shifted by 1,
			// this way ii==0 is recognized correctly
			idxposU[jj]=ii+1; 
			// val= L(jj,:)*D(1:bi,1:bi)^{-1}*L(jj,:)^H
			//    = BL->valE[i_bldl](ii,:)*buff(ii+bi,:)^H
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#else
#ifdef _USE_MKL_
			DISTDOT(&val, &bi, pL,&ldc, pb,&cntL);
#else
			val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#endif
#endif			
			ierr=0;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			sb=buffDgl[jj];
			if (sb<=0.0)
			   ierr=1;
			else
			   rval=sqrt(sb*fabs(val));
#else
			sb=buffDgl[jj].r;
			if (sb<=0.0)
			   ierr=1;
			else
			   rval=sqrt(sb*FABS(val));
#endif
			if (!ierr && rval>0.0) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   // lump diagonal entry
			   if (perturbation>0)
			     buffDgl[jj]+=rval;
			   else {
			     buffDgl[jj]-=rval;
			     if (buffDgl[jj]<=0.0)
			        ierr=1;
			   }
			   // distribute normalization to L(jj,:)*D(1:bi,1:bi)^{-1} 
			   // and L(jj,:) simultaneously
			   val=1.0/sqrt(rval);
#else
			   // lump diagonal entry
			   if (perturbation>0)
			      buffDgl[jj].r+=rval;
			   else {
			      buffDgl[jj].r-=rval;
			      if (buffDgl[jj].r<=0.0)
			         ierr=1;
			   }
			   val.r=1.0/sqrt(rval);val.i=0.0;
#endif
			   // rescale L(jj,:)*D(1:bi,1:bi)^{-1}
			   SCAL(&bi,&val,pL,&ldc);
			   // rescale L(jj,:)
			   SCAL(&bi,&val,pb,&cntL);
#ifdef _PROFILING_
			   time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
			   timeBegin=omp_get_wtime();
#endif
			   // update D(1:bi,1:bi) <- 1*L(jj,:)^H*L(jj,:) + 1*D(1:bi,1:bi)
			   jj=1;
			   GEMM(transa,transb, &bi,&bi,&jj, &alpha,
			        pb,&cntL, pb,&cntL, &beta, 
				AvalD[i_bldl],&bi, 1,1);
#ifdef _PROFILING_
			   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
			   timeBegin=omp_get_wtime();
#endif
			} // end if !ierr and rval>0

			if (ierr) {
			   // re-adjust idxL
			   idxL-=bi;

			   ierr=1;
			   printf("BICHOL encountered zero or negative diagonal part\n");

			   // sort original matrix back such that row indices and their
			   // values are taken in increasing order
			   for (j=0; j<n; j++)
			       QSORT(A->val[j],A->rowind[j],idxL,&(A->ncol[j]));
			   
			   // printf("1.0\n"); fflush(stdout);
			   idxL--;idxL=FREE(idxL);
			   // printf("1.01\n"); fflush(stdout);
			   idxU--;idxU=FREE(idxU);
			   // printf("1.02\n"); fflush(stdout);
			   idxbuff =FREE(idxbuff);
			   // printf("1.03\n"); fflush(stdout);
			   idxposL =FREE(idxposL);
			   // printf("1.04\n"); fflush(stdout);
			   idxposU =FREE(idxposU);
			   // printf("1.05\n"); fflush(stdout);
			   invblock=FREE(invblock);
			   // printf("1.06\n"); fflush(stdout);
			   Ahead   =FREE(Ahead);
			   // printf("1.07\n"); fflush(stdout);
			   Alist   =FREE(Alist);
			   // printf("1.08\n"); fflush(stdout);
			   Afirst  =FREE(Afirst);
			   // printf("1.09\n"); fflush(stdout);
			   BLhead  =FREE(BLhead);
			   // printf("1.010\n"); fflush(stdout);
			   BLlist  =FREE(BLlist);
			   // printf("1.011\n"); fflush(stdout);
			   BLfirst =FREE(BLfirst);
			   // printf("1.012\n"); fflush(stdout);
			   AvalL   =FREE(AvalL);
			   // printf("1.013\n"); fflush(stdout);
			   Arowind =FREE(Arowind);
			   // printf("1.014\n"); fflush(stdout);
			   buffC   =FREE(buffC);
			   // printf("1.015\n"); fflush(stdout);
			   buff    =FREE(buff);
			   // printf("1.016\n"); fflush(stdout);
			   buffD   =FREE(buffD);
			   // printf("1.1\n"); fflush(stdout);
			   for (ii=0; ii<nblocks; ii++) {
			       BL->colind[ii] =FREE(BL->colind[ii] );
			       BL->rowind[ii] =FREE(BL->rowind[ii] );
			       BiD->colind[ii]=FREE(BiD->colind[ii]);
			       BL->valE[ii]   =FREE(BL->valE[ii]   );
			       BiD->valD[ii]  =FREE(BiD->valD[ii]  );
			       AvalD[ii]      =FREE(AvalD[ii]      );
			   }
			   // printf("1.2\n"); fflush(stdout);
			   AvalD   =FREE(AvalD);
			   
			   BL->nblockcol =FREE(BL->nblockcol);
			   BL->nblockrow =FREE(BL->nblockrow);
			   BL->colind    =FREE(BL->colind);
			   BL->rowind    =FREE(BL->rowind);
			   BL->valE      =FREE(BL->valE);
			   
			   BiD->nblockcol=FREE(BiD->nblockcol);
			   BiD->colind   =FREE(BiD->colind   );
			   BiD->valD     =FREE(BiD->valD     );
			   
			   BL->nc=BL->nr=BiD->nc=BiD->nr=0;
			   BL->nblocks=BiD->nblocks=0;
			   // printf("1.3\n"); fflush(stdout);
			   
			   if (flag_ildl1t) {
			      free(blocksize);
			      blocksize=myblocksize;
			      nblocks  =mynblocks;
			   } // end if flag_ildl1t
			   // printf("1.4\n"); fflush(stdout);
	      
			   if (perturbation)
			      buffDgl=FREE(buffDgl);
			   // printf("1.5\n"); fflush(stdout);

			   return (ierr);
			} // end if ierr

		     } // end if-else |pL[k]|>=droptol
		 } // end for ii
		 // re-adjust idxL
		 idxL-=bi;
	      } // end if BL->nblockrow[i_bldl]<ldc


	      // if the "low-rank" update is close to a "full-rank" update,
	      // better recompute the inverse D^{-1} from D directly
	      // similarly, recompute L*D^{-1}
	      if (3*cntU>=bi || !invert_blocks) {
#ifdef _PROFILING_
		 time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
		 timeBegin=omp_get_wtime();
#endif
		 // copy updated AvalD back to BiD->valD[i_bldl]
		 memcpy(BiD->valD[i_bldl],AvalD[i_bldl],(size_t)bi*bi*sizeof(FLOAT));
		 // factorize and invert diagonal block using LAPACK
		 POTRF("l",&bi, BiD->valD[i_bldl],&bi, &ierr,1);
		 if (ierr) {
		    if (ierr<0)
		       printf("LAPACK's POTRF: %ld-th argument had an illegal value\n", -ierr);
		    else
		       printf("LAPACK's POTRF routine encountered zero column in step %ld\n", ierr);

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
		    AvalL   =FREE(AvalL);
		    Arowind =FREE(Arowind);
		    buffC   =FREE(buffC);
		    buff    =FREE(buff);
		    buffD   =FREE(buffD);

		    for (i=0; i<nblocks; i++) {
		        BL->colind[i] =FREE(BL->colind[i] );
			BL->rowind[i] =FREE(BL->rowind[i] );
			BiD->colind[i]=FREE(BiD->colind[i]);
			BL->valE[i]   =FREE(BL->valE[i]   );
			BiD->valD[i]  =FREE(BiD->valD[i]  );
			AvalD[i]      =FREE(AvalD[i]      );
		    }
		    AvalD   =FREE(AvalD);
	      
		    BL->nblockcol =FREE(BL->nblockcol);
		    BL->nblockrow =FREE(BL->nblockrow);
		    BL->colind    =FREE(BL->colind);
		    BL->rowind    =FREE(BL->rowind);
		    BL->valE      =FREE(BL->valE);

		    BiD->nblockcol=FREE(BiD->nblockcol);
		    BiD->colind   =FREE(BiD->colind   );
		    BiD->valD     =FREE(BiD->valD     );
	      
		    BL->nc=BL->nr=BiD->nc=BiD->nr=0;
		    BL->nblocks=BiD->nblocks=0;
	      
		    if (flag_ildl1t) {
		       free(blocksize);
		       blocksize=myblocksize;
		       nblocks  =mynblocks;
		    } // end if flag_ildl1t
	      
		    if (perturbation)
		       buffDgl=FREE(buffDgl);
	      
		    return (ierr);
		 } // end if

		 
		 if (invert_blocks) {
		    POTRI("l",&bi, BiD->valD[i_bldl],&bi, &ierr,1);
		    if (ierr) {
		       if (ierr<0)
			  printf("LAPACK's POTRI: %ld-th argument had an illegal value\n", -ierr);
		       else
			  printf("LAPACK's POTRI routine encountered zero column in step %ld\n", ierr);
		       
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
		       AvalL   =FREE(AvalL);
		       Arowind =FREE(Arowind);
		       buffC   =FREE(buffC);
		       buff    =FREE(buff);
		       buffD   =FREE(buffD);

		       for (i=0; i<nblocks; i++) {
			   BL->colind[i] =FREE(BL->colind[i] );
			   BL->rowind[i] =FREE(BL->rowind[i] );
			   BiD->colind[i]=FREE(BiD->colind[i]);
			   BL->valE[i]   =FREE(BL->valE[i]   );
			   BiD->valD[i]  =FREE(BiD->valD[i]  );
			   AvalD[i]      =FREE(AvalD[i]      );
		       }
		       AvalD   =FREE(AvalD);
	      
		       BL->nblockcol =FREE(BL->nblockcol);
		       BL->nblockrow =FREE(BL->nblockrow);
		       BL->colind    =FREE(BL->colind);
		       BL->rowind    =FREE(BL->rowind);
		       BL->valE      =FREE(BL->valE);

		       BiD->nblockcol=FREE(BiD->nblockcol);
		       BiD->colind   =FREE(BiD->colind   );
		       BiD->valD     =FREE(BiD->valD     );
		    
		       BL->nc=BL->nr=BiD->nc=BiD->nr=0;
		       BL->nblocks=BiD->nblocks=0;
	      
		       if (flag_ildl1t) {
			  free(blocksize);
			  blocksize=myblocksize;
			  nblocks  =mynblocks;
		       } // end if flag_ildl1t
	      
		       if (perturbation)
			 buffDgl=FREE(buffDgl);
	      
		       return (ierr);
		    } // end if

		    // copy strict lower triangular part to strict upper triangular part
		    pC=BiD->valD[i_bldl];
		    pb=pC+1;
		    mm=1;
		    for (ii=1; ii<bi; pb++,ii++) {
		        pD=pC+bi*ii;
			COPY(&ii, pb,&bi, pD,&mm);
			// conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
			for (m=0; m<ii; pD++,m++)
			    pD->i=-pD->i;
#endif
		    } // end for ii
		 } // end if invert_blocks
		 else {
		    // nothing, no inversion, no copying
		 }


		 if (invert_blocks) {
		    // multiply strict lower triangular block by the inverse diagonal block
		    // and copy it to BL
		    // use level-3-BLAS 
		    // C=BL(:,i), A=subdgl(buff)=subdgl([BiD(i,i);BL(:,i)]), B=BiD(i,i)
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
			    BiD->valD[i_bldl],&bi, &beta, BL->valE[i_bldl],&ldc,1,1);
		 } // end if invert_blocks
		 else { // Cholesky case
		    // forward/backward solve using the Cholesky factorization of the diagonal block
		    // copy subdiagonal part of buff to BL since solve works in-place with BL^H
		    m=bi*ldc;
		    // buffC = subdgl(buff)^H = subdgl([BiD(i,i);BL(:,i)])^H
		    pb=buff+bi;
		    pL=buffC;
		    i2=1;
		    for (ii=0; ii<ldc; ii++,pb++,pL+=bi) 
		        COPY(&bi,pb,&cntL,pL,&i2);
		    // conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
		    pL=buffC;
		    for (ii=0; ii<m; ii++,pL++) 
		        pL->i=-pL->i;
#endif
		    // forward/backward solves BiD(i,i)*BiD(i,i)^H * BL^H  <- BL^H
		    //                    <=>  BiD(i,i)*BiD(i,i)^H * buffC <- buffC
		    POTRS("l", &bi, &ldc, BiD->valD[i_bldl],&bi, buffC,&bi, &ierr,1);
		    if (ierr<0) {
		       printf("LAPACK's POTRS: %ld-th argument had an illegal value\n", -ierr);
		       return (ierr);
		    }
		    // BL(:,i)=buffC^H
		    pb=buffC;
		    pL=BL->valE[i_bldl];
		    i2=1;
		    for (ii=0; ii<ldc; ii++,pb+=bi,pL++) 
		        COPY(&bi,pb,&i2,pL,&ldc);
		    // conjugation for the complex case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
		    pL=BL->valE[i_bldl];
		    for (ii=0; ii<m; ii++,pL++) 
		        pL->i=-pL->i;
#endif
		 } // end if-else invert_blocks

		 
		 // downdate the remaining diagonal entries
		 pL=BL->valE[i_bldl];
		 pb=buff+bi;
		 idxL+=bi; // advance idxL temporarily
		 for (ii=0; ii<ldc; ii++,pL++,pb++) {
		     // get true index
		     l=idxL[ii];
		     // a nonzero index indicates a drop index (since the
		     // position is shifted by 1), therefore zero indices
		     // refer to those that are not dropped
		     if (!idxposU[l]) {
		        // val=L(l,:)*D(1:bi,1:bi)^{-1}*L(l,:)^H
		        //    =CONJG(pL[k]^H*pb[k])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		        val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#else
#ifdef _USE_MKL_
		        DISTDOT(&val, &bi, pL,&ldc, pb,&cntL);
#else			
		        val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#endif
#endif
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			buffDgl[l]-=MAX(val,0.0);
#else
			buffDgl[l].r-=MAX(val.r,0.0);
#endif
		     } // end if
		 } // end for ii
		 
		 // compress BL->valE[i_bldl] to those entries that are not dropped
		 // source
		 pb=BL->valE[i_bldl];
		 // destination, source and destination may overlap
		 pL=pb;
		 // compress column-by-column
		 for (jj=0; jj<bi; jj++) {
		     for (ii=0; ii<ldc; ii++,pb++) {
		         // get true index
		         m=idxL[ii];
			 // shift large entries
			 // a nonzero index indicates a drop index (since the
			 // position is shifted by 1), therefore zero indices
			 // refer to those that are not dropped
			 if (!idxposU[m])
			    *pL++=*pb;
		     } // end for ii
		 } // end for jj

		 // (re-)allocate space for indices and lower triangular part
		 l=ldc-cntU;
		 BL->nblockrow[i_bldl]=l;
		 BL->rowind[i_bldl]=(integer *)malloc((size_t)l*sizeof(integer));
		 BL->valE[i_bldl]=(FLOAT *)realloc(BL->valE[i_bldl],(size_t)l*bi*sizeof(FLOAT));
		 
		 // now extract indices
		 pi=BL->rowind[i_bldl];
		 for (ii=0; ii<ldc; ii++) {
		     // get true index
		     jj=idxL[ii];
		     // extract associated indices
		     // a nonzero index indicates a drop index (since the
		     // position is shifted by 1), therefore zero indices
		     // refer to those that are not dropped
		     if (!idxposU[jj]) 
		        *pi++=jj;
		 } // end for ii
		 idxL-=bi; // re-adjust idxL
		 
	      }
	      else { // low-rank update case and invert_blocks
#ifdef _PROFILING_
		 time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
		 timeBegin=omp_get_wtime();
#endif
		 // generate coupling system 
		 // C=I+[L(jj,:)]_jj*D(1:bi,1:bi)^{-1}*[L(jj,:)]_jj^H for 
		 // downdating D(1:bi,1:bi)^{-1} using the Sherman-Morrison-
		 // Woodbury formula and store it in buffC
		 if (cntU*cntU+cntU>nbuffC) {
		    nbuffC=ELBOW_BUFF*(cntU*cntU+cntU)+BUFF_EXT;
		    buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
		 } // end if

		 // init coupling system with I
		 jj=cntU*cntU;
		 pC=buffC;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 for (ii=0; ii<jj; ii++)
		     buffC[ii]=0.0;
		 for (ii=0; ii<cntU; ii++,pC+=cntU+1)
		     *pC=1.0;
#else
		 for (ii=0; ii<jj; ii++)
		     buffC[ii].r=buffC[ii].i=0.0;
		 for (ii=0; ii<cntU; ii++,pC+=cntU+1)
		     pC->r=1.0;
#endif
		 
		 // augment with diagonal and off-diagonal scalar products
		 // rows
		 pL=BL->valE[i_bldl];
		 // columns
		 pb=buff+bi; // remember to skip the diagonal block
		 // target matrix
		 pC=buffC;
		 for (ii=0; ii<cntU; ii++) {
		     // advance to the diagonal position of the coupling system
		     pC+=ii;
		     // get true drop index (column)
		     l=idxU[ii];
		     // where is it physically located inside BL->valE[i_bldl]?
		     k=idxposU[l]-1; // remember shift!
		     // start from the diagonal position
		     for (jj=ii; jj<cntU; jj++) {
		         // get second true drop index (row)
		         l=idxU[jj];
			 // where is it physically located inside BL->valE[i_bldl]?
			 m=idxposU[l]-1; // remember shift
			 // val=L(l_old,:)*D(1:bi,1:bi)^{-1}*L(l_current,:)^H
			 //    =CONJG(pL[k]^H*pb[m])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			 val=DISTDOT(&bi, pL+k,&ldc, pb+m,&cntL);
#else
#ifdef _USE_MKL_
			 DISTDOT(&val, &bi, pL+k,&ldc, pb+m,&cntL);
#else
			 val=DISTDOT(&bi, pL+k,&ldc, pb+m,&cntL);
#endif
#endif
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			 *pC++ +=val;
#else
			 pC->r += val.r;
			 pC->i  =-val.i;
			 pC++;
#endif
		     } // end for jj
		 } // end for ii
		 // factorize coupling system Lc*Lc^H=C
		 if (cntU) {
		    POTRF("l",&cntU, buffC,&cntU, &ierr,1);
		    if (ierr) {
		       if (ierr<0)
			  printf("LAPACK's POTRF(coupling system): %ld-th argument had an illegal value\n", -ierr);
		       else
			  printf("LAPACK's POTRF(coupling system) routine encountered zero column in step %ld\n", ierr);

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
		       AvalL   =FREE(AvalL);
		       Arowind =FREE(Arowind);
		       buffC   =FREE(buffC);
		       buff    =FREE(buff);
		       buffD   =FREE(buffD);
		       
		       for (i=0; i<nblocks; i++) {
			   BL->colind[i] =FREE(BL->colind[i] );
			   BL->rowind[i] =FREE(BL->rowind[i] );
			   BiD->colind[i]=FREE(BiD->colind[i]);
			   BL->valE[i]   =FREE(BL->valE[i]   );
			   BiD->valD[i]  =FREE(BiD->valD[i]  );
			   AvalD[i]      =FREE(AvalD[i]      );
		       }
		       AvalD   =FREE(AvalD);
		       
		       BL->nblockcol =FREE(BL->nblockcol);
		       BL->nblockrow =FREE(BL->nblockrow);
		       BL->colind    =FREE(BL->colind);
		       BL->rowind    =FREE(BL->rowind);
		       BL->valE      =FREE(BL->valE);
		       
		       BiD->nblockcol=FREE(BiD->nblockcol);
		       BiD->colind   =FREE(BiD->colind   );
		       BiD->valD     =FREE(BiD->valD     );
		       
		       BL->nc=BL->nr=BiD->nc=BiD->nr=0;
		       BL->nblocks=BiD->nblocks=0;
	      
		       if (flag_ildl1t) {
			  free(blocksize);
			  blocksize=myblocksize;
			  nblocks  =mynblocks;
		       } // end if flag_ildl1t
		       
		       if (perturbation)
			  buffDgl=FREE(buffDgl);
	      
		       return (ierr);
		    } // end if ierr

		    // solve system Lc^{-1}*[L(l,:)*D(1:bi,1:bi)^{-1}]_l, l drop index 
		    // auxiliary buffer to gather column jj of [L(l,:)*D(1:bi,1:bi)^{-1}]_l
		    pC=buffC+cntU*cntU;
		    // L*D(1:bi,1:bi)^{-1}
		    pL=BL->valE[i_bldl];
		    for (jj=0; jj<bi; jj++,pL+=ldc) {
		        // 1. step gather
		        // gather column jj of [L(l,:)*D(1:bi,1:bi)^{-1}]_l
		        for (ii=0; ii<cntU; ii++) {
			    // get true drop index
			    l=idxU[ii];
			    // where is it physically located inside BL->valE[i_bldl]?
			    k=idxposU[l]-1; // remember shift
			    // move to the associated row inside BL->valE[i_bldl]
			    pC[ii]=pL[k];
			} // end for ii
			// 2. step solve
			ii=1;
			// pC=L_c^{-1}pC
			TRSV("l","n","n", &cntU, buffC,&cntU, pC,&ii, 1,1,1);
			// 3. step scatter
			// scatter column jj back to [L(l,:)*D(1:bi,1:bi)^{-1}]_l
			for (ii=0; ii<cntU; ii++) {
			    // get true drop index
			    l=idxU[ii];
			    // where is it physically located inside BL->valE[i_bldl]?
			    k=idxposU[l]-1; // remember shift
			    // move to the associated row
			    pL[k]=pC[ii];
			} // end for ii
		    } // end for jj

		 
		    // downdate D(1:bi,1:bi)^{-1} using Sherman-Morrison-Woodbury
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    transa="t"; transb="n";
		    alpha=-1.0;
		    beta = 1.0;
#else
		    transa="c"; transb="n";
		    alpha.r=-1.0; alpha.i=0.0;
		    beta.r = 1.0; beta.i =0.0;
#endif
		    pL=BL->valE[i_bldl];
		    for (ii=0; ii<cntU; ii++) {
		        // get true drop index
		        l=idxU[ii];
			// where is it physically located inside BL->valE[i_bldl]?
			k=idxposU[l]-1; // remember shift
			// downdate D(1:bi,1:bi)^{-1} <-
			//          -1*D(1:bi,1:bi)^{-1}[L(l,:)]_l^H*L_C^{-H}
			//            *L_C^{-1}*[L(l,:)]_l*D(1:bi,1:bi)^{-1}
			//          + 1*D(1:bi,1:bi)^{-1}
			// i.e. we have D(1:bi,1:bi)^{-1} <-
			//              -1*BL(l,:)^H*BL(l,:) + 1*D(1:bi,1:bi)^{-1}
			// for all l, l drop indices
			jj=1;
			GEMM(transa,transb, &bi,&bi,&jj, &alpha,
			     pL+k,&ldc, pL+k,&ldc, &beta, 
			     BiD->valD[i_bldl],&bi, 1,1);
		    } // end for ii
#ifdef _PROFILING_
	   time_diagonal_block+=omp_get_wtime()-timeBegin;
#endif

	
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
		    

		 } // end if cntU
		 
		 // downdate those rows of BL which are not dropped since
		 // D(1:bi,1:bi)^{-1} was downdated using Sherman-Morrison-Woodbury
		 pL=BL->valE[i_bldl];
		 pb=buff+bi; // original entries are still in "buff", skip dgl blk
		 pC=pL;
		 idxL+=bi;
		 for (ii=0; ii<ldc; ii++,pL++,pb++) {
		     // get true index
		     jj=idxL[ii];
		     // keep large entries and downdate
		     // a nonzero index indicates a drop index (since the position is shifted
		     // by 1), therefore zero indices refer to those that are not dropped
		     if (!idxposU[jj]) {
		        for (jj=0; jj<cntU; jj++) {
			    // get true drop index
			    l=idxU[jj];
			    // where is it physically located inside BL->valE[i_bldl]?
			    k=idxposU[l]-1; // remember shift
			    // L(jj_old,:)*D(1:bi,1:bi)^{-1} <-
			    // L(jj_old,:)*D(1:bi,1:bi)^{-1} -
			    // L(jj_old,:)*D(1:bi,1:bi)^{-1}*[L(l,:)]_l^H*L_c^{-H}
			    // *L_c^{-1}*[L(l,:)]_l*D(1:bi,1:bi)^{-1}
			    // physically, row vectors:
			    // pL <- pL - pb * pC[k]^H*pC[k]
			    // 1. step val=CONJG(pC[k])*pb^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			    val=DISTDOT(&bi, pC+k,&ldc, pb,&cntL);
#else
#ifdef _USE_MKL_
			    DISTDOT(&val, &bi, pC+k,&ldc, pb,&cntL);
#else
			    val=DISTDOT(&bi, pC+k,&ldc, pb,&cntL);
#endif
#endif
			    // 2. step pL <- (-val)*pC[k] + pL 
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			    val=-val;
#else
			    val.r=-val.r; val.i=-val.i;
#endif
			    AXPY(&bi, &val, pC+k,&ldc, pL,&ldc);
			} // end for jj

			// renew true index
			jj=idxL[ii];
			// also downdate associated diagonal entry 
			// dgl_entry[jj] -=  L(jj,:)*D(1:bi,1:bi)^{-1}*L(jj,:)^H
			//               -=  pb*pL^H
			// 1. step:  val = CONJG(pL)*pb^T
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#else
#ifdef _USE_MKL_
			DISTDOT(&val, &bi, pL,&ldc, pb,&cntL);
#else			
			val=DISTDOT(&bi, pL,&ldc, pb,&cntL);
#endif
#endif
			// 2. step: downdate diagonal entry
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			buffDgl[jj]-=MAX(val,0.0);
#else
			buffDgl[jj].r-=MAX(val.r,0.0);
#endif		     
		     } // if
		 } // end for ii
		 idxL-=bi;

		 // space required for the remaining indices
		 l=ldc-cntU;
		 BL->nblockrow[i_bldl]=l;
		 BL->rowind[i_bldl]=(integer *)malloc((size_t)l*sizeof(integer));
	      
		 // now shift and copy back
		 pL=BL->valE[i_bldl];
		 pb=buff;   // remaining entries are temporarily buffered in "buff"
		 pi=BL->rowind[i_bldl];
		 idxL+=bi;
		 for (ii=0; ii<ldc; ii++,pL++) {
		     // get true index
		     jj=idxL[ii];
		     // shift large entries and extract associated indices
		     // a nonzero index indicates a drop index (since the position is shifted
		     // by 1), therefore zero indices refer to those that are not dropped
		     if (!idxposU[jj]) {
		        *pi++=jj;
			COPY(&bi, pL,&ldc, pb,&l);
			pb++;
		     } // if
		 } // end for ii
		 idxL-=bi;
		 pL=BL->valE[i_bldl];
		 pb=buff;
		 for (ii=0; ii<l; ii++,pL++,pb++) {
		     COPY(&bi, pb,&l, pL,&l);
		 } // for ii
		 BL->valE[i_bldl]=(FLOAT *)realloc(BL->valE[i_bldl],(size_t)l*bi*sizeof(FLOAT));
	      } // end if-else 3*cntU>=bi || !invert_blocks
	      

	      // reset idxposU
	      for (m=0; m<cntU; m++) {
		  t=idxU[m];
		  idxposU[t]=0;
	      } // end for m
	      cntU=0;
	   }
	   else { // no perturbation
	      // apply dropping to the strict lower triangular block of BL
	      pL=BL->valE[i_bldl];
	      pb=buff;   // remaining entries are temporarily buffered in "buff"
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
	      BL->nblockrow[i_bldl]=l;
	      BL->rowind[i_bldl]=(integer *)malloc((size_t)l*sizeof(integer));
	      pL=BL->valE[i_bldl];
	      pb=buff;
	      pi=BL->rowind[i_bldl];
	      for (ii=0; ii<l; ii++,pL++,pb++) {
		  *pi++=idxU[ii];
		  COPY(&bi, pb,&ldc, pL,&l);
	      } // for ii
	      BL->valE[i_bldl]=(FLOAT *)realloc(BL->valE[i_bldl],(size_t)l*bi*sizeof(FLOAT));
	   } // end if-else perturbation
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
	    printf("7.idxposL[%ld]=%ld\n",m,idxposL[m]);
	    fflush(stdout);
	  } // end if
	} // end for m	
#endif
#ifdef _PROFILING_
	time_LDL_update_pass_2+=omp_get_wtime()-timeBegin;
#endif	      
	// END compute block column i of BL using linked list for the rows of BL
	// ---------------------------------------------------------------------
#ifdef PRINT_INFO
	if (invert_blocks)
	   printf("computed inverse diagonal block %ld of BiD\n",i_bldl);
	else
	   printf("factorized diagonal block %ld of BiD\n",i_bldl);
	fflush(stdout);
	printf("        ");
	for (m=0; m<bi; m++)
	    printf("%12ld",BiD->colind[i_bldl][m]);
	printf("\n");fflush(stdout);
	for (ii=0; ii<bi; ii++) {
	    printf("%12ld",BiD->colind[i_bldl][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BiD->valD[i_bldl][m*bi+ii]);
	    }
	    printf("\n");fflush(stdout);
#else
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BiD->valD[i_bldl][m*bi+ii].r);
	    }
	    printf("\n            ");
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BiD->valD[i_bldl][m*bi+ii].i);
	    }
	    printf("\n");fflush(stdout);
#endif
	}
	printf("\n");fflush(stdout);
	printf("computed block column %ld of BL\n",i_bldl);fflush(stdout);
	printf("        ");
	for (m=0; m<bi; m++)
	    printf("%12ld",BiD->colind[i_bldl][m]);
	printf("\n");fflush(stdout);
	for (ii=0; ii<BL->nblockrow[i_bldl]; ii++) {
	    printf("%12ld",BL->rowind[i_bldl][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BL->valE[i_bldl][m*BL->nblockrow[i_bldl]+ii]);
	    }
	    printf("\n");fflush(stdout);
#else
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BL->valE[i_bldl][m*BL->nblockrow[i_bldl]+ii].r);
	    }
	    printf("\n            ");
	    for (m=0; m<bi; m++) {
	        printf("%12.4le",BL->valE[i_bldl][m*BL->nblockrow[i_bldl]+ii].i);
	    }
	    printf("\n");fflush(stdout);
#endif
	}
	printf("\n");fflush(stdout);
#endif

  
	
	


#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// check fill-in of current and previous block in order to find out
	// whether a progressive aggregation is recommended
	if (i_bldl && PROGRESSIVE_AGGREGATION && pa) {
	   // first compute old fill
	   k=BiD->nblockcol[i_bldl-1];
	   l=BL->nblockrow[i_bldl-1];
	   m=BL->nblockrow[i_bldl-1];
	   ii=k*(k+l+m);
	   r=BiD->nblockcol[i_bldl];
	   s=BL->nblockrow[i_bldl];
	   t=BL->nblockrow[i_bldl];
	   ii+=r*(r+s+t);

	   // now determine new fill
	   // checkmark fill of the current sub-diagonal block
	   pi=BL->rowind[i_bldl];
	   for (jj=0; jj<s; jj++) {
	       idxL[cntL]=*pi;
	       idxposL[*pi++]=++cntL;
	   } // end for jj
	   // compute additional fill from the previous sub-diagonal block
	   jj=BLfirst[i_bldl-1];
	   pi=BL->rowind[i_bldl-1]+jj;
	   for (; jj<l; jj++) {
	       // only insert additional fill
	       if (!idxposL[*pi]) {
		  idxL[cntL]=*pi;
		  idxposL[*pi]=++cntL;
	       }
	       pi++;
	   } // end for jj
	   
	   // now we are able to compute the new fill
	   jj=(k+r)*(k+r+2*cntL);



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
	      // printf("block column %2ld: old fill %12ld, new fill %12ld\n",i_bldl,ii-cnt_fill,jj);

	      // incorporate additional fill-in
	      cnt_fill+=jj-ii;
	      // leading dimension of the enlarged BiD->valD[i_bldl]
              ldb=k+r;

	      // enlarge blocks
	      // 1. BiD
	      // BiD will now have some diagonal block of type
	      // [D11_new D12_new]
	      // [D21_new D22    ]
	      // with blocks to be computed
	      BiD->colind[i_bldl]=(integer *)realloc(BiD->colind[i_bldl],
						     (size_t)ldb*sizeof(integer));
	      BiD->valD[i_bldl]  =(FLOAT *)  realloc(BiD->valD[i_bldl],
						     (size_t)ldb*ldb*sizeof(FLOAT));

	      AvalD[i_bldl]      =(FLOAT *)  realloc(AvalD[i_bldl],
						     (size_t)ldb*ldb*sizeof(FLOAT));

	      // 2. BL
	      // BL->colind[i_bldl-1] and BL->valE[i_bldl-1] are interpreted
	      // as some kind of block matrix
	      // [L21]
	      // [L31]
	      // where L21 is associated with block i_bldl (i.e. D22)
	      // and L31 refers to the bottom part
	      // BL->colind[i_bldl] and BL->valE[i_bldl] are now extended to cover
	      // a bigger sub-diagonal block of type
	      // [L31_new  L32]
	      // where L31_new is to be computed from the old block structures
	      // and L32 might be interlaced with a few zero rows
	      BL->colind[i_bldl]=(integer *)realloc(BL->colind[i_bldl],
						    (size_t)ldb*sizeof(integer));
	      BL->rowind[i_bldl]=(integer *)realloc(BL->rowind[i_bldl],
						    (size_t)cntL*sizeof(integer));
	      BL->valE[i_bldl]  =(FLOAT *)  realloc(BL->valE[i_bldl],
						    (size_t)cntL*ldb*sizeof(FLOAT));

	      
	      // we have to keep in mind that at the moment, BL==L, BiD==D
	      // is block unit lower triangular, i.e. we have
	      // a local factorization of type
	      // in the case of invert_blocks:
	      //   [ I   0 ] [D11^{-1}   0     ] [I L21^T L31^T]
	      //   [L21  I ] [  0      D22^{-1}] [0   I   L32^T]
	      //   [L31 L32]
	      // respectively in the other case we have
	      //   [ I   0 ] [LD11   0 ][LD11   0 ]^T [I L21^T L31^T]
	      //   [L21  I ] [  0  LD22][  0  LD22]   [0   I   L32^T]
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
	      pL=BL->valE[i_bldl]+s*r-1;
	      // future position of the last entry of the extended L32 inside
	      // the enlarged block
	      pb=BL->valE[i_bldl]+cntL*ldb-1;
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
	      pL=BL->valE[i_bldl-1]+l*k-1;
	      // future position of the last entry of the extended L31 inside
	      // the enlarged block
	      pb=BL->valE[i_bldl]+cntL*k-1;
	      // index array column i_bldl-1
	      pi=BL->rowind[i_bldl-1];
	      // index array column i_bldl
	      pi2=BL->rowind[i_bldl];
	      // number of rows of L21
	      it=BLfirst[i_bldl-1];
	      for (i2=0; i2<k; i2++,pL-=it) {
	          // position of the last row inside BL->nblockrow[i_bldl-1]
		  ii=l-1; 
	          // position of the last row inside BL->nblockrow[i_bldl]
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
	      pL=BL->valE[i_bldl]+cntL*k;
	      // array of nonzero indices of L21
	      pi=BL->rowind[i_bldl-1];
	      // first index of the diagonal block, we need this information in
	      // order to correctly adjust the pointer to the columns of L32
	      tt=BL->colind[i_bldl][0];
	      // number of columns to be cached
	      mm=BLfirst[i_bldl-1];
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
		      BL->valE[i_bldl-1],&l, &beta,
		      BL->valE[i_bldl],&cntL, 1,1);


	      // 3. finalization
	      // to finalize the aggregated sub-diagonal block, we extend
	      // the column index set as well as the row index set
	      // shift existing indices
	      // old index position last column index
	      pi=BL->colind[i_bldl]+r-1;
	      // new index position last column index
	      pi2=BL->colind[i_bldl]+ldb-1;
	      for (ii=0; ii<r; ii++)
	          *pi2--=*pi--;
	      // insert column indices from the previous block
	      memcpy(BL->colind[i_bldl],BL->colind[i_bldl-1],k*sizeof(integer));
	      BL->nblockcol[i_bldl]=ldb;
	      // insert new row indices in increasing order
	      // position of the last fill row inside idxL
	      ii=cntL-1; 
	      // position of the last existing row inside idxL
	      jj=s-1;
	      // future position of the last index of the extended index set inside
	      // the enlarged block
	      pi=BL->rowind[i_bldl]+cntL-1;
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
	      BL->nblockrow[i_bldl]=cntL;

#ifdef PRINT_INFO
	      printf("extended BL{%ld}\n",i_bldl);fflush(stdout);
	      if (BL->nblockrow[i_bldl]) {
		 printf("        ");
		 for (jj=0; jj<BL->nblockcol[i_bldl]; jj++) {
		     printf("%12ld",BL->colind[i_bldl][jj]);
		 }
		 printf("\n");fflush(stdout);
		 for (ii=0; ii<BL->nblockrow[i_bldl]; ii++) {
		     printf("%12ld",BL->rowind[i_bldl][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		     for (jj=0; jj<BL->nblockcol[i_bldl]; jj++) {
		         printf("%12.4le",BL->valE[i_bldl][BL->nblockrow[i_bldl]*jj+ii]);
		     }
#else
		     for (jj=0; jj<BL->nblockcol[i_bldl]; jj++) {
		         printf("%12.4le",BL->valE[i_bldl][BL->nblockrow[i_bldl]*jj+ii].r);
		     }
		     printf("\n            ");		  
		     for (jj=0; jj<BL->nblockcol[i_bldl]; jj++) {
		         printf("%12.4le",BL->valE[i_bldl][BL->nblockrow[i_bldl]*jj+ii].r);
		     }
#endif
		     printf("\n");	fflush(stdout);	  
		 }
		 printf("\n");fflush(stdout);
	      }
#endif	      
	      // ----------------------------------------------------
	      // END compute aggregated lower triangular block
	      // ----------------------------------------------------

	      

	      // ----------------------------------------------------
	      // compute aggregated regular diagonal block
	      // ----------------------------------------------------	      

	      // [A11_new A12_new]=[  A11        A11*L21^T    ]
	      // [A21_new A22    ] [L21*A11  A22+L21*A11*L21^T]
	      //       ...        =[   A11        A12_new     ]
	      //                   [A12_new^T A22+L21*A12_new ]
	      // 1. shifts & copy

	      // shift current regular diagonal block A22 inside the enlarged block
	      // this will become A22_new^0 which has still to be updated
	      // old position last element
	      pD=AvalD[i_bldl]+r*r-1;
	      // new position last element
	      pb=AvalD[i_bldl]+ldb*ldb-1;
	      // shift and remember different leading dimensions r and ldb
	      for (ii=0; ii<r; ii++,pb-=k)
		  for (jj=0; jj<r; jj++)
		      *pb--=*pD--;
  		      
	      // copy A11 to the leading part of the new enlarged diagonal block
	      // start A11 in the previous block
	      pD=AvalD[i_bldl-1];
	      // start D11_new
	      pb=AvalD[i_bldl];
	      jj=1;
	      for (ii=0; ii<k; ii++,pD+=k,pb+=ldb)
	          // copy column ii of D11
		  COPY(&k, pD,&jj, pb,&jj);

	      // init new super-diagonal block A12_new and sub-diagonal block
	      // A21_new of AvalD with 0
	      // beginning of the super-diagonal block
	      pD=AvalD[i_bldl]+ldb*k;
	      // beginning of the sub-diagonal block
	      pL=AvalD[i_bldl]+k;
	      for (ii=0; ii<r; ii++,pD+=r,pL++) {
		  pU=pL;
		  for (jj=0; jj<k; jj++,pU+=ldb) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
                      *pD++=0.0;
                      *pU=0.0;
#else
		      pD->r=pD->i=0.0;
		      pU->r=pU->i=0.0;
		      pD++;
#endif
                  } // end for jj
	      } // end for ii
	      
              // 2. A12

	      // compute new super-diagonal block using GEMM
	      // A12_new =   A11*CONJG(L21^T)
	      // Since L21 may only have a subset of rows of A21,
	      // we cache the product in buff and scatter the result to
	      // A12_new afterwards
	      // number of rows in BL21
	      ii=BLfirst[i_bldl-1];
	      // use buff to hold the copy
	      if (ii*k>nbuff) {
		 // increase buffer
		 nbuff=ELBOW_BUFF*(ii*k)+BUFF_EXT;
		 buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
	      } // end if
	      // clear target buffer
	      pb=buff;
	      for (jj=0; jj<ii*k; jj++,pb++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  *pb=0.0;
#else
		  pb->r=pb->i=0.0;
#endif
	      } // end for jj
	      // now we are ready to apply GEMM
	      // buff = 1*A11*CONJG(L21^T)+0*buff
	      transa="n"; transb="c";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      alpha=1.0;   beta=0.0;
#else
	      alpha.r=1.0; beta.r=0.0;
	      alpha.i=0.0; beta.i=0.0;
#endif
	      if (k && ii)
		 GEMM(transa,transb, &k,&ii,&k, &alpha,
		      AvalD[i_bldl-1],&k, 
		      BL->valE[i_bldl-1],&m, &beta,
		      buff,&k,1,1);
#ifdef PRINT_INFO3
	      printf("partially extended AvalD{%ld}, before copying\n",i_bldl);fflush(stdout);
	      printf("        ");
	      for (jj=0; jj<ldb; jj++) {
		  printf("%12ld",jj+BiD->colind[i_bldl-1][0]);
	      }
	      printf("\n");fflush(stdout);
	      for (i2=0; i2<ldb; i2++) {
		  printf("%12ld",i2+BiD->colind[i_bldl-1][0]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2]);
		  }
#else
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2].r);
		  }
		  printf("\n            ");
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2].i);
		  }
#endif
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
#endif

	      // scatter buff to A12_new and A21_new
	      // pointer to buff
	      pb=buff;
	      // pointer to A21_new
	      pL=AvalD[i_bldl]+k;
	      // pointer to A12_new
	      pD=AvalD[i_bldl]+ldb*k;
	      // list of nonzero row indices of BL21
	      pi=BL->rowind[i_bldl-1];
	      // first index of the diagonal block, we need this information in
	      // order to correctly adjust the pointer to the columns of A12_new
	      tt=BiD->colind[i_bldl][0];
	      for (jj=0; jj<ii; jj++,pb+=k,pi++) {
		  // scatter column jj of buff to A12_new and A21_new
		  mm=1;
		  COPY(&k, pb,&mm, pD+ldb*(*pi-tt),&mm);
		  COPY(&k, pb,&mm, pL+(*pi-tt),&ldb);
		  // complex conjugation
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_
		  pU=pL+(*pi-tt);
		  for (mm=0; mm<k; mm++,pU+=ldb)
		      pU->i=-pU->i;
#endif
	      } // end for jj
#ifdef PRINT_INFO3
	      printf("partially extended AvalD{%ld}, after copying\n",i_bldl);fflush(stdout);
	      printf("        ");
	      for (jj=0; jj<ldb; jj++) {
		  printf("%12ld",jj+BiD->colind[i_bldl-1][0]);
	      }
	      printf("\n");fflush(stdout);
	      for (i2=0; i2<ldb; i2++) {
		  printf("%12ld",i2+BiD->colind[i_bldl-1][0]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2]);
		  }
#else
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2].r);
		  }
		  printf("\n            ");
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+i2].i);
		  }
#endif
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
#endif

              // 3. A22	      

	      // compute new leading diagonal block
	      // A22_new =   L21*A11*L21^T +   A22
	      //         =   L21*A12_new   +   A22 
	      //         = 1*L21*A12_new   + 1*A22 using GEMM
	      // by pre-multiplying with L21 we only need a small part of
	      // A22 associated with the rows of L21
	      // to do so, copy the required rows and columns of A22 to buffC
	      // "buff" still holds the nonzero columns of A12_new associated
	      // with the nonzero rows of L21
	      // target buffer
	      if (ii*ii>nbuffC) {
		 // increase buffer
		 nbuffC=ELBOW_BUFF*(ii*ii)+BUFF_EXT;
		 buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	      } // end if
	      pb=buffC;
	      // pointer to the array of sub-diagonal row indices in L
	      pi=BL->rowind[i_bldl-1];
	      // pointer to A22
	      pD=AvalD[i_bldl]+ldb*k+k;
	      for (jj=0; jj<ii; jj++) {
		  // column pi[jj]-tt of A22
		  pL=pD+ldb*(pi[jj]-tt);
		  for (i2=0; i2<ii; i2++)
		      *pb++=pL[pi[i2]-tt];
	      } // end for jj
	      
	      // now we are ready to apply GEMM to compute A22_new:
	      // buffC = 1*L21*buff + 1*buffC 
	      transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	      alpha=1.0;   beta=1.0;
#else
	      alpha.r=1.0; beta.r=1.0;
	      alpha.i=0.0; beta.i=0.0;
#endif
	      if (k && ii)
		 GEMM(transa,transb, &ii,&ii,&k, &alpha,
		      BL->valE[i_bldl-1],&l,
		      buff,&k, &beta,
		      buffC,&ii, 1,1);
	      // symmetrize result exactly
	      pb=buffC;
	      for (jj=0; jj<ii; jj++) {
		  // row jj
		  pL=buffC+jj;
#if defined _SINGLE_REAL || defined _DOUBLE_REAL_		  
		  for (i2=0; i2<jj; i2++,pL+=ii)
		      *pL=*pb++;
#else
		  for (i2=0; i2<jj; i2++,pL+=ii) {
		      pL->r= pb->r;
		      pL->i=-pb->i; // conjugate complex value
		      pb++;
		  } // end for i2
		  // diagonal entry is real
		  pL->i=0.0;
#endif
		  pb+=ii-jj;
	      } // end for jj
	      // scatter result back to A22_new
	      pb=buffC;
	      for (jj=0; jj<ii; jj++) {
		  // column pi[jj]-tt of A22
		  pL=pD+ldb*(pi[jj]-tt);
		  for (i2=0; i2<ii; i2++)
		      pL[pi[i2]-tt]=*pb++;
	      } // end for jj
	      
#ifdef PRINT_INFO
	      printf("extended AvalD{%ld}\n",i_bldl);fflush(stdout);
	      printf("        ");
	      for (jj=0; jj<ldb; jj++) {
		  printf("%12ld",jj+BiD->colind[i_bldl-1][0]);
	      }
	      printf("\n");fflush(stdout);
	      for (ii=0; ii<ldb; ii++) {
		  printf("%12ld",ii+BiD->colind[i_bldl-1][0]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+ii]);
		  }
#else
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+ii].r);
		  }
		  printf("\n            ");
		  for (jj=0; jj<ldb; jj++) {
		      printf("%12.4le",AvalD[i_bldl][ldb*jj+ii].i);
		  }
#endif
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
#endif
	      // ----------------------------------------------------
	      // END compute aggregated regular diagonal block
	      // ----------------------------------------------------	      


	      
	      
	      // ----------------------------------------------------
	      // compute aggregated inverse/factorized diagonal block
	      // ----------------------------------------------------	      

	      // inverse case:
	      // [D11_new D12_new]=[D11+L21^T*D22*L21  -L21^T*D22]
	      // [D21_new D22    ] [   -D22*L21            D22   ]
	      //       ...        =[D11-D12_new*L21      -D12_new]
	      //                   [   -D12_new^T          D22   ]
	      // respectively factorized case:
	      // [LD11      *  ]  =[LD11        * ]
	      // [LD21_new LD22]   [L21*LD11  LD22]
	      
	      
	      // 1. shifts & copy
	      
	      // shift current (inverse) diagonal block D22 inside the enlarged block
	      // old position last element
	      pD=BiD->valD[i_bldl]+r*r-1;
	      // recall that in the 1x1 case we do we use 1/D22 rather than sqrt(D22)
	      if (!invert_blocks && r==1) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 *pD=1.0/sqrt(*pD);
#else
		 pD->r=1.0/sqrt(pD->r);
#endif
	      } // end if !invert_blocks and r=1
	      // new position last element
	      pb=BiD->valD[i_bldl]+ldb*ldb-1;
	      // shift and remember different leading dimensions r and ldb
	      for (ii=0; ii<r; ii++,pb-=k)
		  for (jj=0; jj<r; jj++)
		      *pb--=*pD--;
  		      
	      // copy D11 to the leading part of the new enlarged diagonal block
	      // this will become D11_new^0 which has still to be updated
	      // start D11 in the previous block
	      pD=BiD->valD[i_bldl-1];
	      // recall that in the 1x1 case we do we use 1/D11 rather than sqrt(D11)
	      if (!invert_blocks && k==1) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 *pD=1.0/sqrt(*pD);
#else
		 pD->r=1.0/sqrt(pD->r);
#endif
	      } // end if !invert_blocks and k=1
	      // start D11_new
	      pb=BiD->valD[i_bldl];
	      jj=1;
	      for (ii=0; ii<k; ii++,pD+=k,pb+=ldb)
		  // copy column ii of D11
		  COPY(&k, pD,&jj, pb,&jj);

	      // init new super-diagonal block D12_new of BiD with 0
	      // beginning of the super-diagonal block
	      pD=BiD->valD[i_bldl]+ldb*k;
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
		 // D12_new =   -L21^T*D22
		 //         =   -BL21^T*buff
		 // we need to compute D22^T*BL21=D22*BL21 in order to build the
		 // extended inverse diagonal block
		 // extract the rows of D22 that are associated with the nonzero
		 // rows of BL21, copy them to buff, then apply GEMM
		 // number of rows in BL21
		 ii=BLfirst[i_bldl-1];
		 // start of D22 in the extended block
		 pD=BiD->valD[i_bldl]+(ldb+1)*k;
		 // use buff to hold the copy
		 if (ii*r>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*(ii*r)+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 // target buffer
		 pb=buff;
		 // pointer to the array of sub-diagonal row indices
		 pi=BL->rowind[i_bldl-1];
		 // first index of the diagonal block, we need this information in
		 // order to correctly adjust the pointer to the rows of D22
		 tt=BiD->colind[i_bldl][0];
		 for (jj=0; jj<ii; jj++,pb++,pi++)
		     // copy logical row *pi of D22 to buff
		     COPY(&r, pD+(*pi-tt),&ldb, pb,&ii);
	      
/*	      
printf("buff copy:\n");
pb=buff;	      
for (jj=0; jj<ii; jj++,pb++){
    for (i2=0; i2<r; i2++)
        printf("%12.4le",pb[ii*i2]);
    printf("\n");
}
*/
	         // now we are ready to apply GEMM
		 // D12_new =    CONJG(BL21^T)*D22
		 //         = -1*CONJG(BL21^T)*buff + 0*D12_new
		 transa="c"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=-1.0;   beta=0.0;
#else
		 alpha.r=-1.0; beta.r=0.0;
		 alpha.i= 0.0; beta.i=0.0;
#endif
		 if (k && r && ii)
		    GEMM(transa,transb, &k,&r,&ii, &alpha,
			 BL->valE[i_bldl-1],&m, 
			 buff,&ii, &beta,
			 BiD->valD[i_bldl]+ldb*k,&ldb,1,1);
      
/*	      
printf("D11:\n");
pb=BiD->valD[i_bldl-1];
for (jj=0; jj<k; jj++,pb++){
    for (i2=0; i2<k; i2++)
        printf("%12.4le",pb[k*i2]);
    printf("\n");
}
printf("-D11 buffC:\n");
pb=BiD->valD[i_bldl]+ldb*k;	      
for (jj=0; jj<k; jj++,pb++){
    for (i2=0; i2<r; i2++)
        printf("%12.4le",pb[ldb*i2]);
    printf("\n");
}
*/	      

                 // 3. D11	      
		 
		 // compute new leading diagonal block
		 // D11_new =   L21^T*D22*L21 + D11
		 //         =   -D12_new *L21 +   D11 
		 //         = -1*D12_new *L21 + 1*D11 using GEMM
		 // by post-multiplying with L21 we only need a small part of
		 // D12_new associated with the rows of L21
		 // to do so, copy the required columns of D12_new to buff
		 // number of rows in L21
		 ii=BLfirst[i_bldl-1];
		 // start of D12_new
		 pD=BiD->valD[i_bldl]+ldb*k;
		 // target buffer
		 if (k*ii>nbuff) {
		    // increase buffer
		    nbuff=ELBOW_BUFF*(k*ii)+BUFF_EXT;
		    buff=(FLOAT *)realloc(buff,(size_t)nbuff*sizeof(FLOAT));
		 } // end if
		 pb=buff;
		 // pointer to the array of sub-diagonal row indices in L
		 pi=BL->rowind[i_bldl-1];
		 // first index of the diagonal block, we need this information in
		 // order to correctly adjust the pointer to the columns of D12_new
		 tt=BiD->colind[i_bldl][0];
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
			 BL->valE[i_bldl-1],&l, &beta,
			 BiD->valD[i_bldl],&ldb, 1,1);
		 // symmetrize result exactly
		 pb=BiD->valD[i_bldl];
		 for (jj=0; jj<k; jj++) {
		     // row jj
		     pL=BiD->valD[i_bldl]+jj;
#if defined _SINGLE_REAL || defined _DOUBLE_REAL_		  
		     for (i2=0; i2<jj; i2++,pL+=ldb)
		         *pL=*pb++;
#else
		     for (i2=0; i2<jj; i2++,pL+=ldb) {
		         pL->r= pb->r;
			 pL->i=-pb->i; // conjugate complex value
			 pb++;
		     } // end for i2
		     // diagonal entry is real
		     pL->i=0.0;
#endif
		     pb+=ldb-jj;
		 } // end for jj
	      

		 // 4. D21

		 // copy conjugate transposed super-diagonal block to sub-diagonal block 
		 // D21_new =  CONJG(D12_new^T)
		 for (jj=0; jj<r; jj++) {
		     mm=1;
		     pD=BiD->valD[i_bldl]+k+jj;
		     //                                       ->  
		     COPY(&k, BiD->valD[i_bldl]+ldb*(k+jj),&mm, pD,&ldb);
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
		     // conjugate imaginary part
		     for (mm=0; mm<k; pD+=ldb,mm++)
		         pD->i=-pD->i;
#endif		  
		 } // end for jj
	      } // end if invert_blocks
	      else {
		 // compute D21_new=L21*tril(L11) using level-2 BLAS
		 // D21_new(:,mm) = 1*L21(:,mm:end)*D11(mm:end,mm) + 0*D21_new(:,mm)
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
		 pi=BL->rowind[i_bldl-1];
		 // first index of the diagonal block D22, we need this information in
		 // order to correctly adjust the pointer to the rows of D22
		 tt=BiD->colind[i_bldl][0];
		 // start of BL21
		 pD=BL->valE[i_bldl-1];
		 // number of rows in BL21 (=beginning of B31)
		 ii=BLfirst[i_bldl-1];
		 for (jj=0; jj<ii; jj++,pD++,pi++)
		     // copy logical row *pi of BL21 to buff
		     COPY(&k, pD,&l, pb+(*pi-tt),&r);

		 transa="n"; 
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 alpha=1.0;   beta=0.0;
#else
		 alpha.r=1.0; beta.r=0.0;
		 alpha.i=0.0; beta.i=0.0;
#endif
		 jj=1;
		 // start of BD11(old)
		 pD=BiD->valD[i_bldl-1];
		 // start of BD21_new
		 pb=BiD->valD[i_bldl]+k;
		 for (mm=0; mm<k; mm++) {
		     ii=k-mm;
		     GEMV(transa, &r,&ii, &alpha, buff+mm*r,&r,
			  pD+k*mm+mm,&jj, &beta,  pb+mm*ldb,&jj, 1);
		 } // end for mm
	      } // end if-else invert_blocks


	      // 5. finalization
	      
	      // to finalize the aggregated diagonal block, we extend
	      // the column index set
	      // first shift existing indices
	      // old index position last column index
	      pi=BiD->colind[i_bldl]+r-1;
	      // new index position last column index
	      pi2=BiD->colind[i_bldl]+ldb-1;
	      for (ii=0; ii<r; ii++)
		  *pi2--=*pi--;
	      // now insert indices from the previous block
	      memcpy(BiD->colind[i_bldl],BiD->colind[i_bldl-1],k*sizeof(integer));
	      BiD->nblockcol[i_bldl]=ldb;

#ifdef PRINT_INFO
	      printf("extended BiD{%ld}\n",i_bldl);fflush(stdout);
	      printf("        ");
	      for (jj=0; jj<BiD->nblockcol[i_bldl]; jj++) {
		  printf("%12ld",BiD->colind[i_bldl][jj]);
	      }
	      printf("\n");fflush(stdout);
	      for (ii=0; ii<BiD->nblockcol[i_bldl]; ii++) {
		  printf("%12ld",BiD->colind[i_bldl][ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		  for (jj=0; jj<BiD->nblockcol[i_bldl]; jj++) {
		      printf("%12.4le",BiD->valD[i_bldl][BiD->nblockcol[i_bldl]*jj+ii]);
		  }
#else
		  for (jj=0; jj<BiD->nblockcol[i_bldl]; jj++) {
		      printf("%12.4le",BiD->valD[i_bldl][BiD->nblockcol[i_bldl]*jj+ii].r);
		  }
		  printf("\n            ");
		  for (jj=0; jj<BiD->nblockcol[i_bldl]; jj++) {
		      printf("%12.4le",BiD->valD[i_bldl][BiD->nblockcol[i_bldl]*jj+ii].i);
		  }
#endif
		  printf("\n");fflush(stdout);
	      }
	      printf("\n");fflush(stdout);
#endif
	      
	      // ----------------------------------------------------
	      // END compute aggregated inverse diagonal block
	      // ----------------------------------------------------	      






	      



	      

	      // remove old blocks and adapt BLfirst
	      free(BiD->colind[i_bldl-1]);BiD->colind[i_bldl-1]=NULL;
	      free(BiD->valD[i_bldl-1]);  BiD->valD[i_bldl-1]=NULL;
	      BiD->nblockcol[i_bldl-1]=0;

	      free(BL->rowind[i_bldl-1]);BL->rowind[i_bldl-1]=NULL;
	      free(BL->colind[i_bldl-1]);BL->colind[i_bldl-1]=NULL;
	      free(BL->valE[i_bldl-1]);  BL->valE[i_bldl-1]=NULL;
	      BL->nblockcol[i_bldl-1]=0;
	      BL->nblockrow[i_bldl-1]=0;
	      BLfirst[i_bldl-1]=0;

	      free(AvalD[i_bldl-1]);     AvalD[i_bldl-1]=NULL;
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
	       printf("8.idxposL[%ld]=%ld\n",m,idxposL[m]);
	       fflush(stdout);
	     } // end if
	   } // end for m	
	   for (m=0; m<n; m++) {
	     if (idxposU[m]) {
	       printf("8.idxposU[%ld]=%ld\n",m,idxposU[m]);
	       fflush(stdout);
	     } // end if
	   } // end for m	
#endif
	} // end if i_bldl and PROGRESSIVE_AGGREGATION
#ifdef _PROFILING_
	time_progressive_aggregation+=omp_get_wtime()-timeBegin;
#endif
	

	
#ifdef _PROFILING_
	timeBegin=omp_get_wtime();
#endif
	// finally update linked lists
	// recall that BLfirst is already advanced such
	// that its first entry refers to an index beyond block column i
	
	// beginning of the linked list
	it=BLhead[i_bldl];
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
	BLfirst[i_bldl]=0;
	// block column non-empty
	if (0<BL->nblockrow[i_bldl]) {
	   // get row index
	   t=BL->rowind[i_bldl][0];
	   // associated block
	   tt=invblock[t];
	   // add i to the linked list of block column tt
	   BLlist[i_bldl]=BLhead[tt];
	   BLhead[tt]=i_bldl;
	} // end if
#ifdef PRINT_INFO
	printf("linked list BL\n");fflush(stdout);
	for (m=0; m<BL->nblocks; m++)
	    printf("%12ld",BLhead[m]);
	printf("\n");fflush(stdout);
	for (m=0; m<BL->nblocks; m++)
	    printf("%12ld",BLlist[m]);
	printf("\n");fflush(stdout);
	for (m=0; m<BL->nblocks; m++)
	    printf("%12ld",BLfirst[m]);
	printf("\n");
	fflush(stdout);
#endif
#ifdef _PROFILING_
	time_list_update+=omp_get_wtime()-timeBegin;
#endif
	
	// advance to next block column
	i++;
	i_bldl++;
  } // end while j<n
  // **************************************************************************
  // *****                     END main loop                              *****
  // **************************************************************************


#ifdef _PROFILING_
  timeBegin=omp_get_wtime();
#endif
  // post processing

  // shift blocks in order to erase old empty blocks
  i_bldl=0;
  for (i=0; i<BiD->nblocks; i++) {
      // block size of the current block
      bi=BiD->nblockcol[i];

      // if the block was non-empty
      if (bi!=0) {
	 // shift block structures
	 BiD->nblockcol[i_bldl]=BiD->nblockcol[i];
	 BiD->colind[i_bldl]   =BiD->colind[i];
	 BiD->valD[i_bldl]     =BiD->valD[i];

	 BL->nblockcol[i_bldl]=BL->nblockcol[i];
	 BL->nblockrow[i_bldl]=BL->nblockrow[i];
	 BL->colind[i_bldl]   =BL->colind[i];
	 BL->rowind[i_bldl]   =BL->rowind[i];
	 BL->valE[i_bldl]     =BL->valE[i];

	 i_bldl++;
      } // end if
  } // end for i
  BiD->nblocks=BL->nblocks=i_bldl;
  
  // re-arrange sparse block data structures
  BiD->nblockcol=(integer *) realloc(BiD->nblockcol,(size_t)BiD->nblocks*sizeof(integer));
  BiD->colind   =(integer **)realloc(BiD->colind,   (size_t)BiD->nblocks*sizeof(integer *));
  BiD->valD     =(FLOAT **)  realloc(BiD->valD,     (size_t)BiD->nblocks*sizeof(FLOAT *));

  BL->nblockcol=(integer *) realloc(BL->nblockcol,(size_t)BL->nblocks*sizeof(integer));
  BL->nblockrow=(integer *) realloc(BL->nblockrow,(size_t)BL->nblocks*sizeof(integer));
  BL->colind   =(integer **)realloc(BL->colind,   (size_t)BL->nblocks*sizeof(integer *));
  BL->rowind   =(integer **)realloc(BL->rowind,   (size_t)BL->nblocks*sizeof(integer *));
  BL->valE     =(FLOAT **)  realloc(BL->valE,     (size_t)BL->nblocks*sizeof(FLOAT *));

  
  
  
#ifdef PRINT_INFO
  printf("Final factorization\n");
  printf("BiD, BL\n");fflush(stdout);
  for (i=0; i<BiD->nblocks; i++) {
      printf("block %3ld of BL, BiD\n",i);fflush(stdout);
      printf("diagonal block BiD[%3ld]\n            ",i);
      for (m=0; m<BiD->nblockcol[i]; m++) 
	  printf("%12ld",BiD->colind[i][m]);
      printf("\n");fflush(stdout);
      for (k=0; k<BiD->nblockcol[i]; k++) {
	  printf("%12ld",BiD->colind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%12.4le",BiD->valD[i][m*BiD->nblockcol[i]+k]);
	  printf("\n");fflush(stdout);
#else	    
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%12.4le",BiD->valD[i][m*BiD->nblockcol[i]+k].r);
	  printf("\n            ");
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%12.4le",BiD->valD[i][m*BiD->nblockcol[i]+k].i);
	  printf("\n");fflush(stdout);
#endif
      }
      printf("\n");fflush(stdout);

      printf("sub-diagonal block BL[%3ld]\n            ",i);
      if (BL->nblockrow[i]>0) {
	 for (m=0; m<BL->nblockcol[i]; m++) 
	     printf("%12ld",BL->colind[i][m]);
	 printf("\n");fflush(stdout);
	 for (k=0; k<BL->nblockrow[i]; k++) {
	     printf("%12ld",BL->rowind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%12.4le",BL->valE[i][m*BL->nblockrow[i]+k]);
#else
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%12.4le",BL->valE[i][m*BL->nblockrow[i]+k].r);
	     printf("\n            ");
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%12.4le",BL->valE[i][m*BL->nblockrow[i]+k].i);
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
  free(idxbuff);
  idxL--;free(idxL);
  free(idxposL);
  idxU--;free(idxU);
  free(idxposU);
  free(invblock);
  free(Ahead);
  free(Alist);
  free(Afirst);
  free(BLhead);
  free(BLlist);
  free(BLfirst);

  for (i=0; i<nblocks; i++)
      if (AvalD[i]!=NULL)
	 free(AvalD[i]);
  free(AvalD);
  free(AvalL);
  free(Arowind);
  free(buffD);
  free(buffC);
  free(buff);

  // no a priori block partitioning given, use ILDL(1,droptol)
  if (flag_ildl1t) {
     free(blocksize);
     // restore user's values
     blocksize=myblocksize;
     nblocks  =mynblocks;
  } // end if flag_ildl1t

  if (perturbation)
     buffDgl=FREE(buffDgl);

  
  // finally include diagonal scaling to the determinant
  if (SL!=NULL)
     for (i=0; i<n; i++)
         determinant->r-=2*log(SL[i]);
  // make sure to consider the main branch
  while (determinant->i>M_PI)
        determinant->i-=2*M_PI;
  while (determinant->i<=-M_PI)
        determinant->i+=2*M_PI;
  
  BL->nnz=0;
  for (i=0; i<BL->nblocks-1; i++) 
      BL->nnz+=BL->nblockrow[i]*BL->nblockcol[i];
  BiD->nnz=0;
  for (i=0; i<BiD->nblocks; i++) 
      BiD->nnz+=BiD->nblockcol[i]*BiD->nblockcol[i];
  
#ifdef _PROFILING_
  time_post_processing+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
  time_bilu=omp_get_wtime()-time_bilu;
  printf("profiling summary\n");
  printf("0) initialization matrix and lists            %8.1le\n",time_init_matrix);
  if (ildl1t==BLOCK_ILU1T)
     printf("1) ILDL(1,%8.1le)                           %8.1le\n",droptol,time_ilu1t);
  else if (ildl1t==BLOCK_ILUPT)
     printf("1) ILDL(%2ld,%8.1le)                          %8.1le\n",level_of_fill,droptol,time_ilu1t);
  else if (ildl1t==BLOCK_SUPERNODES)
     printf("1) SUPERNODES                                 %8.1le\n",time_ilu1t);
  else if (ildl1t==BLOCK_APP_SUPERNODES)
     printf("1) approximate SUPERNODES                     %8.1le\n",time_ilu1t);
  else
     printf("1) ---                                        %8.1le\n",time_ilu1t);
  
  printf("2) Extraction of blocks from A                %8.1le\n",time_extract_block);
  printf("3) LDL update first pass                      %8.1le\n",time_LDL_update_pass_1);
  printf("   Scalar LDL update                          %8.1le\n",time_scalar_update);
  printf("   Rank-1 LDL update                          %8.1le\n",time_rank_1_update);
  printf("   Small size block LDL update                %8.1le\n",time_small_block_update);
  printf("   Level-3 BLAS LDL update                    %8.1le\n",time_level_3_blas_update);
  printf("   Linked list update                         %8.1le\n",time_list_update);
  printf("   Remainder LDL update pass 2                %8.1le\n",time_LDL_update_pass_2);
  printf("4) Diagonal block factorization and inversion %8.1le\n",time_diagonal_block);
  printf("5) Off-diagonal block re-scaling and dropping %8.1le\n",time_off_diagonal_block);
  printf("6) Progressive aggregation                    %8.1le\n",time_progressive_aggregation);
  printf("7) Post processing                            %8.1le\n",time_post_processing);
  printf("Total BICHOL time %8.1le\n\n",time_bilu);

  fflush(stdout);
#endif
  
  return (ierr);
} // end bilu
