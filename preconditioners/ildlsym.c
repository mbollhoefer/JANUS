/* $Id: ildlsym.c 4002 2018-02-14 16:17:43Z bolle $ */
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

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif




#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYILDLSYM     SSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       ssytrf
#define MYSYMTRI       ssytri
#elif defined _DOUBLE_REAL_
#define MYILDLSYM     DSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       dsytrf
#define MYSYMTRI       dsytri
#elif defined _SINGLE_COMPLEX_
#define MYILDLSYM     CSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       csytrf
#define MYSYMTRI       csytri
#else // double complex
#define MYILDLSYM     ZSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       zsytrf
#define MYSYMTRI       zsytri
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYILDLSYM     SSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       ssytrf
#define MYSYMTRI       ssytri
#elif defined _DOUBLE_REAL_
#define MYILDLSYM     DSYMilu
#define CONJG(A)       (A)
#define MYSYMTRF       dsytrf
#define MYSYMTRI       dsytri
#elif defined _SINGLE_COMPLEX_
#define MYILDLSYM     CHERilu
#define CONJG(A)       (-(A))
#define MYSYMTRF       chetrf
#define MYSYMTRI       chetri
#else // double complex
#define MYILDLSYM     ZHERilu
#define CONJG(A)       (-(A))
#define MYSYMTRF       zhetrf
#define MYSYMTRI       zhetri
#endif //-if-elif-else single-real

#endif //-else complex-symmetric



#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

// relaxation parameters when re-allocating buffers
#define ELBOW_BUFF      1.2
#define BUFF_EXT        100




// #define PRINT_CHECK 
// #define PRINT_INFO 
// #define printf mexPrintf
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif




integer MYILDLSYM(SPARSEMATRIX A,
		  SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD,
		  SPARSEBLOCKMATRIX *BUT,
		  REALS *SL, REALS *SR, integer *p, integer *invq,
		  doublecomplex *determinant,
		  integer *isdefinite,
		  REALS droptol)
{  
  /*
    compute (block) incomplete LDL^T factorization with relative drop tolerance

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

    Output
    ------
    BL          lower block triangular matrix L with unit diagonal
    BiD         block diagonal matrix D^{-1} (LAPACK-factorized and inverted)
    BUT         not referenced
  */

  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,mm,tt,it,it_next,jt,jjt,jt_next,i2, // counters
          *invblock,   // mapping scalar index j -> block i of A
          nrAL,ncAL,nrncAL, // estimates for the blocks of A
          *pi,*pi2,    // index array pointers
          nbuff,nbuffC,nbuffD,// logical size of buff, buffC, buffD
          lda,ldb,ldc, // leading dimension used in combination with level-3-BLAS
          flag,        // flag indicating when a block column is exceeded
          bi,          // current block size
          nrowind,     // length of array Arowind
          i_bldl,      // index counter for the block ILDL
          cnt_fill,    // counter for additional fill when merging blocks
          *rowind,ncol,*colind,nrow,// short cuts
    
          *idxL,*idxU, // temporary index list of row indices
          *idxposL,*idxposU,// associated inverse mapping storing the position
          cntL,        // length of idxL
          cntL_old,    // old length of idxL when needed

          cntU,        // length of idxU
          cntU_old,    // old length of idxU when needed

          *idxbuff,    // auxiliary index buffer
          skip_gthr_sctr, // flag when gather/scatter operation is not needed

          *Ahead,      // head of the linked list for the rows of A
          *Alist,      // (incomplete) linked list for the rows of A
          *Afirst,     // physical position of the first nonzero entry of each 
                       // column as far as we still need to access this column
                       // to extract the associated rows
    
          *BLhead,     // head of the linked list of block columns of BL
                       // required to extract the scalar rows of BL
          *BLlist,     // (incomplete) linked list of block columns of BL
                       // required to extract the scalar rows of BL
          *BLfirst,    // physical position of the first nonzero entry in each
                       // block column of BL as far as it is required to access
                       // a scalar row of BL
    
          *Arowind,    // array with the nonzero row indices of the current
                       // sub-diagonal block of A

          nblocks,     // number of blocks and and block size
          *blocksize,
    
          ierr=0;      // return value ilu

  
  FLOAT   **AvalD,     // buffer for the diagonal blocks of A
          *AvalL,      // buffer for the current sub-diagonal block of A
          *AvalU,      // buffer for the current transposed super-diagonal block of A
          *pL,*pD,*pU,*pb,*pC, // pointers to lower/diagonal/upper/buff/buffC
          val,valD,    // temporary scalar numerical value
          *buff,       // buffer for the current block column/row of BL
          *buffC,      // buffer for level-3-BLAS update
          *buffD,      // buffer for diagonal block x upper triangular block
          alpha, beta; // parameters for level-3-BLAS

  REALS   rval;        // temporary scalar real numerical value

  char    *transa, *transb; // strings for level-3-BLAS

#ifdef _PROFILING_
  double timeBegin,
         time_ilu,
         time_init_matrix,
         time_extract_block=0.0,
         time_LDL_update_pass_1=0.0,
         time_LDL_update_pass_2=0.0,
         time_scalar_update=0.0,
         time_rank_1_update=0.0,
         time_small_block_update=0.0,
         time_list_update=0.0,
         time_diagonal_block=0.0,
         time_off_diagonal_block=0.0,
         time_progressive_aggregation=0.0,
         time_post_processing;

  time_ilu=omp_get_wtime();
#endif
  
  
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

  
#ifdef PRINT_CHECK
  for (j=0; j<n; j++) {
      if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n)
	 printf("ILDL ERROR, wrong p/invq, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j], invq[p[j]]);
  }  
  if (SL!=NULL) {
     for (j=0; j<n; j++) {
       if (fabs(SL[j])==0.0)
	  printf("ILDL ERROR, zero scaling, step %ld\n", j);
     }  
  }
  for (j=0; j<n; j++) {
      for (k=0; k<A->ncol[j]; k++) {
          i=A->rowind[j][k];
	  if (i<0 || i>=n)
	     printf("ILDL ERROR, illegal matrix index (%ld,%ld)\n", i,j);
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
 
  // set up arrays for A(p,p)^T (A is stored by by columns)
  // implicitly these structures work with A(p,p)
  // head of the linked list for the rows of A(p,p)
  Ahead =(integer *)malloc((size_t)n*sizeof(integer));
  // linked list for the leading rows of A(p,p) up to step j
  Alist =(integer *)calloc((size_t)n,sizeof(integer));
  // position of the first entry at step j
  Afirst=(integer *)calloc((size_t)n,sizeof(integer));
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


  
  
#ifdef _PROFILING_
  timeBegin=omp_get_wtime();
#endif
  // initially only 1x1 pivots
  nblocks=n;
  blocksize=(integer *)malloc((size_t)nblocks*sizeof(integer));
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
  for (i=0; i<nblocks; i++)
      BL->nblockrow[i]=0;

  // set up arrays for the rows BL^T (BL is stored by block columns)
  // these are required to access BL by rows during the computation
  // of the current block column of BL
  // head of the linked list for the rows of BL
  BLhead =(integer *)malloc((size_t)nblocks*sizeof(integer));
  // linked list for the leading rows of BL up to step j
  BLlist =(integer *)calloc((size_t)nblocks,sizeof(integer));
  // position of the first entry at step j
  BLfirst=(integer *)calloc((size_t)nblocks,sizeof(integer));
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
		     pD[jt-r].i=CONJG(val.i);
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
		     pL[tt-1].i=CONJG(val.i);
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
		    val.r=-valD.r*pU->r       +valD.i*CONJG(pU->i);
		    val.i=-valD.r*CONJG(pU->i)-valD.i*pU->r;
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
			   val.r=-valD.r*pU->r       +valD.i*CONJG(pU->i);
			   val.i=-valD.r*CONJG(pU->i)-valD.i*pU->r;
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
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   for (i2=0; i2<l; i2++) {
			       for (mm=0; mm<l; mm++)
				   printf("%20.12le",AvalD[it][l*mm+i2]);
			       printf("\n");fflush(stdout);
			   }
#else
			   for (i2=0; i2<l; i2++) {
			       for (mm=0; mm<l; mm++)
				   printf("%20.12le",AvalD[it][l*mm+i2].r);
			       printf("\n");fflush(stdout);
			   }
			   printf("\n");fflush(stdout);
			   for (i2=0; i2<l; i2++) {
			       for (mm=0; mm<l; mm++)
				   printf("%20.12le",AvalD[it][l*mm+i2].i);
			       printf("\n");fflush(stdout);
			   }
#endif
			   printf("\n");fflush(stdout);
			   printf("BL{%ld}(%ld,:)^T\n",it,idxU[m]);fflush(stdout);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   for (mm=0; mm<l; mm++)
			       printf("%12.4le\n",pU[mm*ldb]);
#else
			   for (mm=0; mm<l; mm++)
			       printf("%12.4le\n",pU[mm*ldb].r);
			   printf("\n");fflush(stdout);
			   for (mm=0; mm<l; mm++)
			       printf("%12.4le\n",pU[mm*ldb].i);
#endif
			   printf("\n");fflush(stdout);
				   
#endif
			   mm=1;
			   // buffD = 1*D_{it,it}*BL{it}(m,:)^T+0*buffD
			   transa="t";
			   GEMV(transa,&l,&l,&alpha, AvalD[it],&l, pU,&ldb,
				&beta, buffD,&mm,1);
#ifdef PRINT_INFO1
			   printf("diagonal block D(%ld,%ld)*BL{%ld}(%ld,:)^T\n",it,it,it,idxU[m]);fflush(stdout);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   for (mm=0; mm<l; mm++) {
			       printf("%20.12le\n",buffD[mm]);
			   }
#else
			   for (mm=0; mm<l; mm++) {
			       printf("%20.12le\n",buffD[mm].r);
			   }
			   printf("\n");fflush(stdout);
			   for (mm=0; mm<l; mm++) {
			       printf("%20.12le\n",buffD[mm].i);
			   }
#endif
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
#ifdef _COMPLEX_SYMMETRIC_
			       val=DISTDOTU(&l, buffD,&mm, pL,&lda);
#else
			       val=DISTDOT(&l, buffD,&mm, pL,&lda);
#endif
			       pb[i2].r-=val.r;
			       pb[i2].i-=val.i;
#endif
			       pL++;
			   } // end for k
		       } // end for m


#ifdef PRINT_INFO
		       printf("temporarily computed block column %ld of BL+D\n",i_bldl);fflush(stdout);
		       printf("        ");
		       for (m=r; m<s; m++)
			   printf("%20ld",m);
		       printf("\n");fflush(stdout);
		       for (ii=0; ii<cntL; ii++) {
			   printf("%20ld",idxL[ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			   for (m=0; m<bi; m++) {
			       printf("%20.12le",buff[m*cntL+ii]);
			   }
			   printf("\n");fflush(stdout);
#else
			   for (m=0; m<bi; m++) {
			       printf("%20.12le",buff[m*cntL+ii].r);
			   }
			   printf("\n                    ");fflush(stdout);
			   for (m=0; m<bi; m++) {
			     printf("%20.12le",buff[m*cntL+ii].i);
			   }
			   printf("\n");fflush(stdout);
#endif
		       }
		       printf("\n");fflush(stdout);
#endif


		       
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
		    transa="n"; transb="c";		    
#else
		    for (mm=0; mm<l*cntU; mm++) 	        
		        buffD[mm].r=buffD[mm].i=0.0;
		    alpha.r= 1.0; alpha.i=0.0;
		    beta.r = 0.0; beta.i =0.0;
#ifdef _COMPLEX_SYMMETRIC_
		    transa="n"; transb="t";
#else
		    transa="n"; transb="c";
#endif
#endif
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

#ifdef PRINT_INFO
		    printf("temporarily computed block column %ld of BL+D\n",i_bldl);fflush(stdout);
		    printf("        ");
		    for (m=r; m<s; m++)
		        printf("%20ld",m);
		    printf("\n");fflush(stdout);
		    for (ii=0; ii<cntL; ii++) {
		        printf("%20ld",idxL[ii]);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
			for (m=0; m<bi; m++) {
			    printf("%20.12le",buff[m*cntL+ii]);
			}
			printf("\n");fflush(stdout);
#else
			for (m=0; m<bi; m++) {
			    printf("%20.12le",buff[m*cntL+ii].r);
			}
			printf("\n                    ");fflush(stdout);
			for (m=0; m<bi; m++) {
			    printf("%20.12le",buff[m*cntL+ii].i);
			}
			printf("\n");fflush(stdout);
#endif
		    }
		    printf("\n");fflush(stdout);
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

	// force diagonal entries to be real-valued in the complex Hermitian case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_ && !defined _COMPLEX_SYMMETRIC_
	pb=buff;
	for (m=0; m<bi; m++,pb+=cntL+1)
	    pb->i=0.0;
#endif
	// for technical reasons, make a copy of the diagonal block
	// (stored in the begining of buff) to AvalD[i_bldl] before inversion
	
	pb=buff; pD=AvalD[i_bldl];
	for (m=0; m<bi; m++,pb+=cntL,pD+=bi) {
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
	   if (rval==0.0) {
	      ierr=1;
	      printf("ILDL encountered zero pivot\n");

	      free(blocksize);
	      // ... more free();
	      return (ierr);
	   }
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
	   rval=sqrt(rval)*droptol; // |diagonal_entry|*droptol
	   pb=buff+1; // start strict lower triangular part
	   l=0;
	   for (ii=0; ii<ldc; ii++,pb++)
	       if (FABS(*pb)>=rval)
		  l++;
	   // now we know the fill
	   BL->nblockrow[i_bldl]=l;
	   BL->valE[i_bldl]=(FLOAT *)malloc((size_t)l*sizeof(FLOAT));
	   BL->rowind[i_bldl]=(integer *)malloc((size_t)l*sizeof(integer));
	   	   
	   // multiply strict lower triangular part by the inverse diagonal entry
	   // and copy it to BL
	   pb=buff+1;
	   pL=BL->valE[i_bldl];
	   pi=BL->rowind[i_bldl];
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
	else {
#ifdef _PROFILING_
	   timeBegin=omp_get_wtime();
#endif
	   ii=1; // optimal block size for SYTRF
	   transa="l";
	   jj=-1;
#if defined _SINGLE_REAL_
	   m=ILAENV(&ii, "SSYTRF", transa, &bi, &jj, &jj, &jj,6,1);
#elif defined _DOUBLE_REAL_
	   m=ILAENV(&ii, "DSYTRF", transa, &bi, &jj, &jj, &jj,6,1);
#elif defined _SINGLE_COMPLEX
#ifdef _COMPLEX_SYMMETRIC_
	   m=ILAENV(&ii, "CSYTRF", transa, &bi, &jj, &jj, &jj,6,1);
#else
	   m=ILAENV(&ii, "CSHERF", transa, &bi, &jj, &jj, &jj,6,1);
#endif
#else
#ifdef _COMPLEX_SYMMETRIC_
	   m=ILAENV(&ii, "ZSYTRF", transa, &bi, &jj, &jj, &jj,6,1);
#else
	   m=ILAENV(&ii, "ZHETRF", transa, &bi, &jj, &jj, &jj,6,1);
#endif
#endif
	   // use buffC as work array
	   l=m*bi;
	   if (l>nbuffC) {
	      nbuffC=ELBOW_BUFF*l+BUFF_EXT;
	      buffC=(FLOAT *)realloc(buffC,(size_t)nbuffC*sizeof(FLOAT));
	   } // end if
	   // factorize and invert diagonal block using LAPACK
	   // use idxU as buffer for pivot indices
	   MYSYMTRF(transa,&bi, buff,&cntL, idxbuff, buffC,&l, &ierr,1);
	   if (ierr) {
	      if (ierr<0)
		 printf("LAPACK's MYSYMTRF: %ld-th argument had an illegal value\n", -ierr);
	      else
		 printf("LAPACK's MYSYMTRF routine encountered zero column in step %ld\n", ierr);

	      free(blocksize);
	      // more free();...
	      return (ierr);
	   } // end if
	   MYSYMTRI("l",&bi, buff,&cntL, idxbuff, buffC, &ierr,1);
	   if (ierr) {
	      if (ierr<0)
		 printf("LAPACK's MYSYMTRI: %ld-th argument had an illegal value\n", -ierr);
	      else
		 printf("LAPACK's MYSYMTRI routine encountered zero column in step %ld\n", ierr);

	      free(blocksize);
	      // more free();...
	      return (ierr);
	   } // end if
	   // copy strict lower triangular part to strict upper triangular part
	   pb=buff+1;
	   mm=1;
	   for (ii=1; ii<bi; pb++,ii++) {
	       pD=buff+cntL*ii;
	       COPY(&ii, pb,&cntL, pD,&mm);
	       // cojugation for the complex Hermitian case
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_ && !defined _COMPLEX_SYMMETRIC_
	       for (m=0; m<ii; pD++,m++)
		   pD->i=-pD->i;
#endif
	   } // end for ii
	   
	   // extract inverse diagonal block to BiD
	   BiD->valD[i_bldl]=(FLOAT *)malloc((size_t)bi*bi*sizeof(FLOAT));
	   pD=BiD->valD[i_bldl];
	   pb=buff;
	   for (ii=0; ii<bi; ii++,pD+=bi,pb+=cntL)
	       memcpy(pD,pb,(size_t)bi*sizeof(FLOAT));
	   
#ifdef PRINT_INFO0
	   printf("computed inverse diagonal block %ld of BiD\n",i_bldl);
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
	   // multiply strict lower triangular block by the inverse diagonal block
	   // and copy it to BL
	   ldc=cntL-bi; // number of rows before dropping
	   k=ldc*bi;    // total size
	   BL->valE[i_bldl]=(FLOAT *)malloc((size_t)k*sizeof(FLOAT));
	   // use level-3-BLAS 
	   // init with 0
	   pL=BL->valE[i_bldl];
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
		   BiD->valD[i_bldl],&bi, &beta, BL->valE[i_bldl],&ldc,1,1);
	
	   // apply dropping to the strict lower triangular block of BL
	   pL=BL->valE[i_bldl];
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
#ifdef _PROFILING_
	   time_off_diagonal_block+=omp_get_wtime()-timeBegin;
#endif
	} // end if-else bi=1

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
	printf("computed inverse diagonal block %ld of BiD\n",i_bldl);fflush(stdout);
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
      printf("diagonal block BiD[%3ld]\n                    ",i);
      for (m=0; m<BiD->nblockcol[i]; m++) 
	  printf("%20ld",BiD->colind[i][m]);
      printf("\n");fflush(stdout);
      for (k=0; k<BiD->nblockcol[i]; k++) {
	  printf("%20ld",BiD->colind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%20.13le",BiD->valD[i][m*BiD->nblockcol[i]+k]);
	  printf("\n");fflush(stdout);
#else	    
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%20.13le",BiD->valD[i][m*BiD->nblockcol[i]+k].r);
	  printf("\n                    ");
	  for (m=0; m<BiD->nblockcol[i]; m++) 
	      printf("%20.13le",BiD->valD[i][m*BiD->nblockcol[i]+k].i);
	  printf("\n");fflush(stdout);
#endif
      }
      printf("\n");fflush(stdout);

      printf("sub-diagonal block BL[%3ld]\n                    ",i);
      if (BL->nblockrow[i]>0) {
	 for (m=0; m<BL->nblockcol[i]; m++) 
	     printf("%20ld",BL->colind[i][m]);
	 printf("\n");fflush(stdout);
	 for (k=0; k<BL->nblockrow[i]; k++) {
	     printf("%20ld",BL->rowind[i][k]);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%20.13le",BL->valE[i][m*BL->nblockrow[i]+k]);
#else
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%20.13le",BL->valE[i][m*BL->nblockrow[i]+k].r);
	     printf("\n                    ");
	     for (m=0; m<BL->nblockcol[i]; m++) 
	         printf("%20.13le",BL->valE[i][m*BL->nblockrow[i]+k].i);
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
  free(blocksize);
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
  free(blocksize);
  // more free(); ...
#ifdef _PROFILING_
  time_post_processing+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
  time_ilu=omp_get_wtime()-time_ilu;
  printf("profiling summary\n");
  printf("0) initialization matrix and lists            %12.4le\n",time_init_matrix);
  printf("1) Extraction of blocks from A                %12.4le\n",time_extract_block);
  printf("2) LDL update first pass                      %12.4le\n",time_LDL_update_pass_1);
  printf("   Scalar LDL update                          %12.4le\n",time_scalar_update);
  printf("   Rank-1 LDL update                          %12.4le\n",time_rank_1_update);
  printf("   Small size block LDL update                %12.4le\n",time_small_block_update);
  printf("   Linked list update                         %12.4le\n",time_list_update);
  printf("   Remainder LDL update pass 2                %12.4le\n",time_LDL_update_pass_2);
  printf("3) Diagonal block factorization and inversion %12.4le\n",time_diagonal_block);
  printf("4) Off-diagonal block re-scaling and dropping %12.4le\n",time_off_diagonal_block);
  printf("5) Post processing                            %12.4le\n",time_post_processing);
  printf("Total ILDL time %12.4le\n\n",time_ilu);

  fflush(stdout);
#endif
  
  return (ierr);
} // end symilu

 

 
