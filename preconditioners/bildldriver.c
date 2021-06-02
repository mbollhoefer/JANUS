/* $Id: bildldriver.c 7315 2021-05-28 21:00:20Z bolle $ */
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

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define TAU      0.8

/* #define PRINT_CHECK  */
/* #define PRINT_INFO   */
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif


#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYSYMBILU      SSYMbilu
#define MYBILDLDRIVER  SSYMbilu_driver
#define MYSYMAMGINIT   SSYMAMGinit
#define MYSYMPERMMC64NULL        SSYMperm_mc64_null
#define MYSYMPERMMC64AMD         SSYMperm_mc64_amd
#define MYSYMPERMMC64METISN      SSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE      SSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS     SSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM         SSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    SSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     SSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  SSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  SSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS SSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     SSYMperm_matching_rcm
#define MYSYMPERMMNULL           SSYMperm_null
#define MYSYMPERMMAMD            SSYMperm_amd
#define MYSYMPERMMMETISN         SSYMperm_metis_n
#define MYSYMPERMMMETISE         SSYMperm_metis_e
#define MYSYMPERMMMTMETIS        SSYMperm_mtmetis
#define MYSYMPERMMRCM            SSYMperm_rcm
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYSYMBILU      DSYMbilu
#define MYBILDLDRIVER  DSYMbilu_driver
#define MYSYMAMGINIT   DSYMAMGinit
#define MYSYMPERMMC64NULL    DSYMperm_mc64_null
#define MYSYMPERMMC64AMD     DSYMperm_mc64_amd
#define MYSYMPERMMC64METISN  DSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE  DSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS DSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     DSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    DSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     DSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  DSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  DSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS DSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     DSYMperm_matching_rcm
#define MYSYMPERMMNULL       DSYMperm_null
#define MYSYMPERMMAMD        DSYMperm_amd
#define MYSYMPERMMMETISN     DSYMperm_metis_n
#define MYSYMPERMMMETISE     DSYMperm_metis_e
#define MYSYMPERMMMTMETIS    DSYMperm_mtmetis
#define MYSYMPERMMRCM        DSYMperm_rcm
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYSYMBILU      CSYMbilu
#define MYBILDLDRIVER  CSYMbilu_driver
#define MYSYMAMGINIT   CSYMAMGinit
#define MYSYMPERMMC64NULL    CSYMperm_mc64_null
#define MYSYMPERMMC64AMD     CSYMperm_mc64_amd
#define MYSYMPERMMC64METISN  CSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE  CSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS CSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     CSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    CSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     CSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  CSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  CSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS CSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     CSYMperm_matching_rcm
#define MYSYMPERMMNULL       CSYMperm_null
#define MYSYMPERMMAMD        CSYMperm_amd
#define MYSYMPERMMMETISN     CSYMperm_metis_n
#define MYSYMPERMMMETISE     CSYMperm_metis_e
#define MYSYMPERMMMTMETIS    CSYMperm_mtmetis
#define MYSYMPERMMRCM        CSYMperm_rcm
#define CONJG(A)       (A)
#else // double complex
#define MYSYMBILU      ZSYMbilu
#define MYBILDLDRIVER  ZSYMbilu_driver
#define MYSYMAMGINIT   ZSYMAMGinit
#define MYSYMPERMMC64NULL    ZSYMperm_mc64_null
#define MYSYMPERMMC64AMD     ZSYMperm_mc64_amd
#define MYSYMPERMMC64METISN  ZSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE  ZSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS ZSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     ZSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    ZSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     ZSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  ZSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  ZSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS ZSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     ZSYMperm_matching_rcm
#define MYSYMPERMMNULL       ZSYMperm_null
#define MYSYMPERMMAMD        ZSYMperm_amd
#define MYSYMPERMMMETISN     ZSYMperm_metis_n
#define MYSYMPERMMMETISE     ZSYMperm_metis_e
#define MYSYMPERMMMTMETIS    ZSYMperm_mtmetis
#define MYSYMPERMMRCM        ZSYMperm_rcm
#define CONJG(A)       (A)
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYSYMBILU      SSYMbilu
#define MYBILDLDRIVER  SSYMbilu_driver
#define MYSYMAMGINIT   SSYMAMGinit
#define MYSYMPERMMC64NULL    SSYMperm_mc64_null
#define MYSYMPERMMC64AMD     SSYMperm_mc64_amd
#define MYSYMPERMMC64METISN  SSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE  SSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS SSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     SSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    SSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     SSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  SSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  SSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS SSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     SSYMperm_matching_rcm
#define MYSYMPERMMNULL       SSYMperm_null
#define MYSYMPERMMAMD        SSYMperm_amd
#define MYSYMPERMMMETISN     SSYMperm_metis_n
#define MYSYMPERMMMETISE     SSYMperm_metis_e
#define MYSYMPERMMMTMETIS    SSYMperm_mtmetis
#define MYSYMPERMMRCM        SSYMperm_rcm
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYSYMBILU      DSYMbilu
#define MYBILDLDRIVER  DSYMbilu_driver
#define MYSYMAMGINIT   DSYMAMGinit
#define MYSYMPERMMC64NULL    DSYMperm_mc64_null
#define MYSYMPERMMC64AMD     DSYMperm_mc64_amd
#define MYSYMPERMMC64METISN  DSYMperm_mc64_metis_n
#define MYSYMPERMMC64METISE  DSYMperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS DSYMperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     DSYMperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    DSYMperm_matching_null
#define MYSYMPERMMATCHINGAMD     DSYMperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  DSYMperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  DSYMperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS DSYMperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     DSYMperm_matching_rcm
#define MYSYMPERMMNULL       DSYMperm_null
#define MYSYMPERMMAMD        DSYMperm_amd
#define MYSYMPERMMMETISN     DSYMperm_metis_n
#define MYSYMPERMMMETISE     DSYMperm_metis_e
#define MYSYMPERMMMTMETIS    DSYMperm_mtmetis
#define MYSYMPERMMRCM        DSYMperm_rcm
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYSYMBILU      CHERbilu
#define MYBILDLDRIVER  CHERbilu_driver
#define MYSYMAMGINIT   CHERAMGinit
#define MYSYMPERMMC64NULL    CHERperm_mc64_null
#define MYSYMPERMMC64AMD     CHERperm_mc64_amd
#define MYSYMPERMMC64METISN  CHERperm_mc64_metis_n
#define MYSYMPERMMC64METISE  CHERperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS CHERperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     CHERperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    CHERperm_matching_null
#define MYSYMPERMMATCHINGAMD     CHERperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  CHERperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  CHERperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS CHERperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     CHERperm_matching_rcm
#define MYSYMPERMMNULL       CHERperm_null
#define MYSYMPERMMAMD        CHERperm_amd
#define MYSYMPERMMMETISN     CHERperm_metis_n
#define MYSYMPERMMMETISE     CHERperm_metis_e
#define MYSYMPERMMMTMETIS    CHERperm_mtmetis
#define MYSYMPERMMRCM        CHERperm_rcm
#define CONJG(A)       (-(A))
#else // double complex
#define MYSYMBILU      ZHERbilu
#define MYBILDLDRIVER  ZHERbilu_driver
#define MYSYMAMGINIT   ZHERAMGinit
#define MYSYMPERMMC64NULL    ZHERperm_mc64_null
#define MYSYMPERMMC64AMD     ZHERperm_mc64_amd
#define MYSYMPERMMC64METISN  ZHERperm_mc64_metis_n
#define MYSYMPERMMC64METISE  ZHERperm_mc64_metis_e
#define MYSYMPERMMC64MTMETIS ZHERperm_mc64_mtmetis
#define MYSYMPERMMC64RCM     ZHERperm_mc64_rcm
#define MYSYMPERMMATCHINGNULL    ZHERperm_matching_null
#define MYSYMPERMMATCHINGAMD     ZHERperm_matching_amd
#define MYSYMPERMMATCHINGMETISN  ZHERperm_matching_metis_n
#define MYSYMPERMMATCHINGMETISE  ZHERperm_matching_metis_e
#define MYSYMPERMMATCHINGMTMETIS ZHERperm_matching_mtmetis
#define MYSYMPERMMATCHINGRCM     ZHERperm_matching_rcm
#define MYSYMPERMMNULL       ZHERperm_null
#define MYSYMPERMMAMD        ZHERperm_amd
#define MYSYMPERMMMETISN     ZHERperm_metis_n
#define MYSYMPERMMMETISE     ZHERperm_metis_e
#define MYSYMPERMMMTMETIS    ZHERperm_mtmetis
#define MYSYMPERMMRCM        ZHERperm_rcm
#define CONJG(A)       (-(A))
#endif //-if-elif-else single-real

#endif //-else complex-symmetric

integer MYBILDLDRIVER(SPARSEMATRIX *C,
		      SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT, integer *pivots,
		      REALS *SL, REALS *SR, integer *p, integer *invq, 
		      doublecomplex *determinant,
		      integer *isdefinite,
		      JanusOptions options)
{  
    integer            i,j,k,l,ii,ll,jj,m,n,*invq2,*p2,*invq3,*p3,r,s,mynblocks,*myblocksize=NULL,
                       nB=0, ierr,bi, *invblock=NULL, *blockind, nnz, 
                       cnt, *idx, *idxpos;
    FLOAT              *pcolscale2, *prowscale2, *colscale, *rowscale;
    size_t             mrows, ncols;
    ILUPACKPARAM       ilupackoptions;
    CSRMAT             A,B;
    SPARSEMATRIX       CC;
#ifdef PRINT_CHECK
    integer *check_buff;
#endif
#ifdef _PROFILING_
    double timeBegin,
           time_total=0.0,
           time_bilu=0.0,
           time_pre_processing=0.0,
           time_cosine=0.0,
           time_reordering=0.0;

    timeBegin=time_total=omp_get_wtime();
#endif

#ifdef PRINT_INFO
    printf("entering MYBILDLDRIVER\n"); fflush(stdout);
#endif
    // recompute nnz
    n=C->nr;
    nnz=0;
    for (j=0; j<n; j++)
        nnz+=C->ncol[j];
#ifdef PRINT_CHECK
    check_buff=(integer *)calloc(n,sizeof(integer));;
#endif
    
    /* copy input matrix to sparse row format */
    A.nc=A.nr=n;
    A.ia=(integer *)malloc((size_t)(n+1)*sizeof(integer));
    A.ja=(integer *)malloc((size_t)nnz  *sizeof(integer));
    A. a=(FLOAT *)  malloc((size_t)nnz  *sizeof(FLOAT));
    
    A.ia[0]=1;
    l=0;
    for (j=0; j<n; j++) {
	for (k=0; k<C->ncol[j]; k++) {
	    i=C->rowind[j][k];
	    A.ja[l]=i+1;
	    A.a [l++]=C->val[j][k];
	}
        A.ia[j+1]=A.ia[j]+C->ncol[j];
    }
    A.nnz=nnz;

    /* initialize ILUPACK options structure to its default options */
    MYSYMAMGINIT(&A,&ilupackoptions);
#ifdef _PROFILING_
    time_pre_processing=omp_get_wtime()-timeBegin;
#endif

    
#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
    /* row scaling vector */
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    colscale=rowscale=SL;
#else
    colscale=(FLOAT *)malloc((size_t)n*sizeof(FLOAT));
    rowscale=colscale;
#endif
    nB=n;
#ifdef _PROFILING_
    time_pre_processing+=omp_get_wtime()-timeBegin;
#endif



    
    
    if (options.nblocks==n||options.nblocks==0||options.blocksize==NULL) {
       /* if matching is requested */
       if (options.matching && options.cosine) {
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* call maximum weight matching followed by a symmetrizing routine */
#ifdef PRINT_INFO
	  printf("SYMBILU: apply matching\n");fflush(stdout);
#endif
	  ierr=MYSYMPERMMC64NULL(A, rowscale,colscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: matching performed (%ld 1x1, %ld 2x2) w.r.t. n=%ld\n",nB,(n-nB)/2,n);fflush(stdout);
#endif
	  if (ierr) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
#else
	     free(colscale);
#endif
	     myblocksize=FREE(myblocksize);
	     free(A.ia);
	     free(A.ja);
	     free(A.a);
	     ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	     ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	     ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	     ilupackoptions.daux =FREE(ilupackoptions.daux);
	     return ierr;
	  } // end if ierr

	  /* now the permuted system A(p,p) refers to a matrix, where the leading
	     nB x nB block refers to 1x1 blocks whereas the remaining (n-nB) x (n-nB) 
	     prefers to work with 2x2 pivots. This block is to be maintained in the
	     sequel.
	     Note also that A has been rescaled by some D_rowscale*A*D_colscale
	  */
	  /* the scalings are real-valued anyway */
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	  for (i=0; i<n; i++) {
              SL[i]=colscale[i].r;
	  } /* end for */
#endif
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  } /* end for i */
#ifdef PRINT_CHECK
	  for (i=0; i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;	      
	  for (i=0; i<n; i++) {
	      if (invq[p[i]]!=i) {
	         printf("SYMBILU ERROR, inverse permutation does not match, step %ld, p[%ld]=%ld, invq[p[%ld]]=%ld\n",i,i,p[i],i,invq[p[i]]);
		 fflush(stdout);
	      }
	  }
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
       
	  /* if we want to use the cosine blocking strategy on top of matching,
	     we first have to replace A by a compressed counter part
	  */

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif

	  /* start by copying the permuted pattern */
	  CC.nc=CC.nr=n;
	  CC.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
	  CC.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
	  CC.val   =(FLOAT **)  malloc((size_t)n*sizeof(FLOAT *));
	     
	  for (j=0; j<n; j++) {
	      ii=p[j];
	      CC.ncol[j]=C->ncol[ii];
	      CC.rowind[j]=(integer *)malloc((size_t)CC.ncol[j]*sizeof(integer));
	      CC.val[j]   =(FLOAT *)  malloc((size_t)CC.ncol[j]*sizeof(FLOAT));
	      for (k=0; k<CC.ncol[j]; k++) {
		  CC.rowind[j][k]=invq[C->rowind[ii][k]];
	      } // end for k
	  } /* end for j */
#ifdef PRINT_INFO
	  /* DPrintMatrix(C); */
	  printf("SYMBILU: permuted sparse matrix data structures set up\n");fflush(stdout);
#endif


	  /* compress matrix, use invq as buffer temporarily */
	  /* the leading nB rows correspond to 1x1 pivots */
	  for (i=0; i<nB; i++)
	      invq[i]=i;
	  j=nB;
	  for (; i<n; i+=2)
	      invq[i]=invq[i+1]=j++;

	  /* overwrite row indices by their compressed representation, possibly leading to
	     duplicate indices 
	  */
	  for (j=0; j<n; j++) {
	      for (k=0; k<CC.ncol[j]; k++) {
		  i=CC.rowind[j][k];
		  CC.rowind[j][k]=invq[i];
	      } /* end for k */
	  } /* end for j */

	  /* use invq temporarily as buffer for flags */
	  for (i=0; i<n; i++)
	      invq[i]=0;

	  /* skip duplicate entries in the leading nB rows */
	  for (j=0; j<nB; j++) {
	      /* new start of column j */
	      k=0;
	      for (i=0; i<CC.ncol[j]; i++) {
		  ii=CC.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     CC.rowind[j][k++]=ii;
		  } /* end if */
	      } /* end for i */
	      /* clear check mark array */
	      for (i=0; i<k; i++) {
		  ii=CC.rowind[j][i];
		  invq[ii]=0;
	      }
	      /* update number of nonzero row indices in column j */
	      CC.ncol[j]=k;
	      /* reduce amount of memory */
	      CC.rowind[j]=(integer *)realloc(CC.rowind[j],(size_t)k*sizeof(integer));
	      CC.val[j]   =(FLOAT *)  realloc(CC.val[j],   (size_t)k*sizeof(FLOAT));
	  } /* end for j */
	  
	  /* merge duplicate entries and rows in the remaining n-nB rows */
	  for (m=nB; j<n; m++,j++) {
	      /* new start of column m */
	      k=0;

	      /* make sure that the current column is long enough */
	      CC.rowind[m]=(integer *)realloc(CC.rowind[m],(size_t)(CC.ncol[j]+CC.ncol[j+1])*sizeof(integer));
	      /* first scan column j */
	      for (i=0; i<CC.ncol[j]; i++) {
		  ii=CC.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     CC.rowind[m][k++]=ii;
		  } /* end if */
	      } /* end for i */

	      j++;
	      /* second scan column j+1 */
	      for (i=0; i<CC.ncol[j]; i++) {
		  ii=CC.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     CC.rowind[m][k++]=ii;
		  } /* end if */
	      } /* end for i */
	     
	      /* clear check mark array of column j and j+1 */
	      for (i=0; i<k; i++) {
		  ii=CC.rowind[m][i];
		  invq[ii]=0;
	      }
	      /* update number of nonzero row indices in column j */
	      CC.ncol[m]=k;
	      /* reduce amount of memory */
	      CC.rowind[m]=(integer *)realloc(CC.rowind[m],(size_t)k*sizeof(integer));
	      CC.val[m]   =(FLOAT *)  realloc(CC.val[m],   (size_t)k*sizeof(FLOAT));
	  } /* end for m */
	  CC.nr=CC.nc=m;

	  /* get rid of the remaining columns */
	  for (j=m; j<n; j++) {
	      free(CC.rowind[j]);
	      free(CC.val[j]);
	  }
#ifdef PRINT_INFO
	  printf("SYMBILU: companion matrix for cosine-based algorithm computed, size m=%ld\n",m);fflush(stdout);
#endif

	  /* apply the cosine-based blocking to the compressed matrix */
	  myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
	  p2         =(integer *)malloc((size_t)n*sizeof(integer));
	  mynblocks=0;
	  for (j=0; j<m; j++)
	      p2[j]=invq[j]=j;
	  MCOSINE_SBLOCKS(&CC, p2,invq, myblocksize, &mynblocks, TAU);
#ifdef PRINT_INFO
	  printf("SYMBILU: cosine applied\n");fflush(stdout);
#endif
#ifdef PRINT_CHECK
	  for (i=0; i<m; i++) {
	      if (p2[i]<0 || p2[i]>=m) {
		 printf("SYMBILU ERROR, wrong permutation, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=m) {
		 printf("SYMBILU ERROR, wrong permutation, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;
	  for (i=0; i<m; i++) {
	      if (check_buff[p2[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
		 fflush(stdout);
	      }
	      check_buff[p2[i]]=-1;
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;
	  for (i=0; i<m; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;	      
	  for (i=0; i<m; i++) {
	      if (invq[p2[i]]!=i) {
	         printf("SYMBILU ERROR, inverse permutation does not match, step %ld, p2[%ld]=%ld, invq[p2[%ld]]=%ld\n",i,i,p2[i],i,invq[p2[i]]);
		 fflush(stdout);
	      }
	  }
	  if (mynblocks<=0 || mynblocks>m) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS: block size mismatch %ld\n",mynblocks);
	     fflush(stdout);
	  }
	  j=0;
	  for (i=0; i<mynblocks; i++) {
	      if (myblocksize[i]<=0 || myblocksize[i]>m) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS: block %ld has an illegal value %ld\n",i,myblocksize[i]);
		 fflush(stdout);
	      }
	      j+=myblocksize[i];
	      if (j>m) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS: step %ld, accumulated block size %ld exceeds limit %ld\n",i,j,m);
		 fflush(stdout);
	      }
	  }
	  if (j!=m) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS: accumulated block size %ld does not match %ld\n",j,m);
	     fflush(stdout);
	  }
#endif      
	  

	  /* prolongate permutation to its original size, use invq as auxiliary vector */
	  // compute inverse mapping index->block
	  invblock=(integer *)malloc((size_t)n*sizeof(integer));
	  j=0;
	  for (i=0; i<mynblocks; i++) {
	      bi=myblocksize[i];
	      r=j; s=r+bi;
	      while (j<s) {
		    invblock[j]=i;
		    j++;
	      } /* end while */
	  } /* end for i */

	  // expand blocks of B(p2,p2)
	  j=0;
	  for (i=0;i<m; i++) {
	      // consider B(:,k), B(k,:)
	      k=p2[i];
	      /* 1x1 block, no expansion */
	      if (k<nB)
		 invq[j++]=k;
	      else { /* 2x2 block, k refers to a block of two rows/columns */
		 myblocksize[invblock[i]]++;
		 invq[j]=2*k-nB;
		 invq[j+1]=invq[j]+1;
		 j+=2;
	      }
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong expanded permutation, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  if (mynblocks<=0 || mynblocks>n) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS expanded: block size mismatch %ld\n",mynblocks);
	     fflush(stdout);
	  }
	  j=0;
	  for (i=0; i<mynblocks; i++) {
	      if (myblocksize[i]<=0 || myblocksize[i]>n) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS expanded: block %ld has an illegal value %ld\n",i,myblocksize[i]);
		 fflush(stdout);
	      }
	      j+=myblocksize[i];
	      if (j>n) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS expanded: step %ld, accumulated block size %ld exceeds limit %ld\n",i,j,n);
		 fflush(stdout);
	      }
	  }
	  if (j!=n) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS expanded: accumulated block size %ld does not match %ld\n",j,n);
	     fflush(stdout);
	  }
#endif      
	  

	  /* build product permutation for the columns (symmetric maximum 
	     weighted matching followed cosine-based blocking)
	  */
	  for (i=0; i<n; i++) 
	      p2[i]=p[invq[i]];
	  /* rewrite permutation */
	  for (i=0; i<n; i++)
	      p[i]=p2[i];
	  /* associated inverse permutation */
	  for (i=0; i<n; i++)
	      invq[p[i]]=i;
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  /* remove remaining columns */
	  for (j=0; j<m; j++) {
	      free(CC.rowind[j]);
	      free(CC.val[j]);
	  }
	  free(p2);
	  free(CC.ncol);
	  free(CC.rowind);
	  free(CC.val);
	  
#ifdef PRINT_INFO
	  printf("SYMBILU: permutation after cosine expanded to incorporate 2x2 pivots\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_cosine+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
#ifdef PRINT_INFO
	  printf("SYMBILU: proper sparse matrix data structures based on A set up\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* now we have applied maximum weight matching followed by cosine-based blocking,
	     yet, a symmetric reordering taking into account the induced block structure
	     still has to be applied
	  */
       } /* end if matching and cosine */
       else if (options.matching) { /* only matching */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  // if (!strcmp("mtmetis",options.ordering)) {
	  if (options.ordering==PERM_MTMETIS) {
#ifdef PRINT_INFO
	     printf("SYMBILU: apply matching+mtmetis\n");fflush(stdout);
#endif
#ifdef _MC64_MATCHING_
	     ierr=MYSYMPERMMC64MTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	     ierr=MYSYMPERMMATCHINGMTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("SYMBILU: matching+mtmetis performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  else if (options.ordering==PERM_METISN){ /* METIS nested dissection by nodes */
#ifdef PRINT_INFO
	     printf("SYMBILU: apply matching+metisn\n");fflush(stdout);
#endif
#ifdef _MC64_MATCHING_
	     ierr=MYSYMPERMMC64METISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	     ierr=MYSYMPERMMATCHINGMETISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif	     
#ifdef PRINT_INFO
	     printf("SYMBILU: matching+metisn performed\n");fflush(stdout);
#endif
	  }
	  //else if (!strcmp("metise",options.ordering)) {
	  else if (options.ordering==PERM_METISE) {
#ifdef PRINT_INFO
	     printf("SYMBILU: apply matching+metise\n");fflush(stdout);
#endif
#ifdef _MC64_MATCHING_
	     ierr=MYSYMPERMMC64METISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	     ierr=MYSYMPERMMATCHINGMETISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("SYMBILU: matching+metise performed\n");fflush(stdout);
#endif
	  }
#endif // _USE_METIS4_
	  // else if (!strcmp("amd",options.ordering)) {
	  else if (options.ordering==PERM_AMD) {
#ifdef _MC64_MATCHING_
	     ierr=MYSYMPERMMC64AMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
             ierr=MYSYMPERMMATCHINGAMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("SYMBILU: matching+amd performed\n");fflush(stdout);
#endif
	  }
	  // else (!strcmp("rcm",options.ordering)) {
	  else {
#ifdef _MC64_MATCHING_
	     ierr=MYSYMPERMMC64RCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	     ierr=MYSYMPERMMATCHINGRCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("SYMBILU: matching+rcm performed\n");fflush(stdout);
#endif
	  }
	  if (ierr) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
#else
	     free(colscale);
#endif
	     myblocksize=FREE(myblocksize);
	     free(A.ia);
	     free(A.ja);
	     free(A.a);
	     ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	     ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	     ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	     ilupackoptions.daux =FREE(ilupackoptions.daux);
	     return ierr;
	  } // end if ierr

#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
	      SL[i]=colscale[i].r;
	  } /* end for */
#endif
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  /* no cosine blocking */
	  myblocksize=NULL;
	  mynblocks=0;

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* we are ready to apply BILDL */
       } /* end if only matching */
       else if (options.cosine) { /* only cosine */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do some scaling in advance */
	  ierr=MYSYMPERMMNULL(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
	  if (ierr) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
#else
	     free(colscale);
#endif
	     free(A.ia);
	     free(A.ja);
	     free(A.a);
	     ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	     ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	     ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	     ilupackoptions.daux =FREE(ilupackoptions.daux);
	     return ierr;
	  } // end if ierr
	  /* the scalings are real-valued anyway */
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	  for (i=0; i<n; i++) {
	      SL[i]=colscale[i].r;
	  } /* end for */
#endif
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif

	  /* apply the cosine-based blocking to the compressed matrix */
	  myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
	  mynblocks=0;

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  MCOSINE_SBLOCKS(C, p,invq, myblocksize,&mynblocks, TAU);
#ifdef PRINT_INFO
	  printf("SYMBILU: cosine applied\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_cosine+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  if (mynblocks<=0 || mynblocks>n) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS: block size mismatch %ld\n",mynblocks);
	     fflush(stdout);
	  }
	  j=0;
	  for (i=0; i<mynblocks; i++) {
	      if (myblocksize[i]<=0 || myblocksize[i]>n) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS: block %ld has an illegal value %ld\n",i,myblocksize[i]);
		 fflush(stdout);
	      }
	      j+=myblocksize[i];
	      if (j>n) {
		 printf("SYMBILU ERROR, MCOSINE_SBLOCKS: step %ld, accumulated block size %ld exceeds limit %ld\n",i,j,n);
		 fflush(stdout);
	      }
	  }
	  if (j!=n) {
	     printf("SYMBILU ERROR, MCOSINE_SBLOCKS: accumulated block size %ld does not match %ld\n",j,n);
	     fflush(stdout);
	  }
#endif      
	  /* we have applied the cosine blocking only and still need to reorder
	     the system taking into account the associated block structure 
	  */
       } /* end if only cosine */
       else { /* neither cosine nor matching nor blocks */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do scaling and simple reordering without asking for matchings */
	  // if (!strcmp("mtmetis",options.ordering)) {
	  if (options.ordering==PERM_MTMETIS) {
	     ierr=MYSYMPERMMMTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: mtmetis performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  // else if (!strcmp("metisn",options.ordering)) {
	  else if (options.ordering==PERM_METISN) {
	     ierr=MYSYMPERMMMETISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: metisn performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("metise",options.ordering)) {
	  else if (options.ordering==PERM_METISE) {
	     ierr=MYSYMPERMMMETISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: metise performed\n");fflush(stdout);
#endif
	  }
#endif // _USE_METIS4_
	  // else if (!strcmp("amd",options.ordering)) {
	  else if (options.ordering==PERM_AMD) {
	     ierr=MYSYMPERMMAMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: amd performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("rcm",options.ordering)) {
	  else if (options.ordering==PERM_RCM) {
	     ierr=MYSYMPERMMRCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: rcm performed\n");fflush(stdout);
#endif
	  }
	  else { /* none */
	     ierr=MYSYMPERMMNULL(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("SYMBILU: no symmetric reordering used\n");fflush(stdout);
#endif
	  }
	  /* the scalings are real-valued anyway */
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	  for (i=0; i<n; i++) {
	      SL[i]=colscale[i].r;
	  } /* end for */
#endif
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif

	  /* no cosine blocking */
	  myblocksize=NULL;
	  mynblocks=0;
	  
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* we are ready to apply BILDL */
       }
    } /* end if no given blocks */
    else { /* block partitioning prescribed */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do some scaling in advance */
	  ierr=MYSYMPERMMNULL(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
	  /* the scalings are real-valued anyway */
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	  for (i=0; i<n; i++) {
	      SL[i]=colscale[i].r;
	  } /* end for */
#endif
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 printf("SYMBILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  myblocksize=options.blocksize;
	  mynblocks=options.nblocks;
#ifdef PRINT_CHECK
	  if (mynblocks<=0 || mynblocks>n) {
	     printf("SYMBILU ERROR, block size mismatch %ld\n",mynblocks);
	     fflush(stdout);
	  }
	  j=0;
	  for (i=0; i<mynblocks; i++) {
	      if (myblocksize[i]<=0 || myblocksize[i]>n) {
		 printf("SYMBILU ERROR, block %ld has an illegal value %ld\n",i,myblocksize[i]);
		 fflush(stdout);
	      }
	      j+=myblocksize[i];
	      if (j>n) {
		 printf("SYMBILU ERROR, step %ld, accumulated block size %ld exceeds limit %ld\n",i,j,n);
		 fflush(stdout);
	      }
	  }
	  if (j!=n) {
	     printf("SYMBILU ERROR, accumulated block size %ld does not match %ld\n",j,n);
	     fflush(stdout);
	  }
#endif      

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* before starting BILDL, we are going to reorder the system based on the
	     given block structure 
	  */
    } /* end if-else no given blocks */

    /* matrix A is not longer required */
    free(A.ia);  
    free(A.ja);
    free(A.a);


    /* if only maximum weight matching was requested without cosine and if no
       initial block partitioning was passed (similar case for 
       no matching+no cosine+no blocks), then we are already done 
       otherwise (cosine or given blocks) we need to apply the symmetric 
       reordering to some compressed companion matrix
    */
    if (myblocksize!=NULL) {
#ifdef PRINT_CHECK
       if (mynblocks<=0 || mynblocks>n) {
	  printf("SYMBILU ERROR, block size mismatch %ld\n",mynblocks);
	  fflush(stdout);
       }
       j=0;
       for (i=0; i<mynblocks; i++) {
	   if (myblocksize[i]<=0 || myblocksize[i]>n) {
	      printf("SYMBILU ERROR, block %ld has an illegal value %ld\n",i,myblocksize[i]);
	      fflush(stdout);
	   }
	   j+=myblocksize[i];
	   if (j>n) {
	      printf("SYMBILU ERROR, step %ld, accumulated block size %ld exceeds limit %ld\n",i,j,n);
	      fflush(stdout);
	   }
       }
       if (j!=n) {
	  printf("SYMBILU ERROR, accumulated block size %ld does not match %ld\n",j,n);
	  fflush(stdout);
       }
#endif      
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* compute a companion matrix B of A(p,p) compressing the associated blocks */
       if (invblock==NULL)
	  invblock=(integer *)malloc((size_t)n*sizeof(integer));
       j=0;
       for (i=0; i<mynblocks; i++) {
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         invblock[j]=i;
		 j++;
	   } /* end while */
       } /* end for i */

       B.nr=B.nc=mynblocks;
       B.ia  =(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
       idx   =(integer *)malloc((size_t)mynblocks*sizeof(integer));
       idxpos=(integer *)calloc((size_t)mynblocks,sizeof(integer));
       cnt=0;
       
       /* first pass: detect nnz column by column */
       j=0;
       nnz=0;
       for (i=0; i<mynblocks; i++) {
	   /* scan block column i */
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         jj=p[j]; /* shortcut */
		 for (k=0; k<C->ncol[jj]; k++) {
		     ll=C->rowind[jj][k];
		     /* we will have q[l]=ll, so in total we are considering C(q[l],p[j]) */
		     l=invq[ll];
		     m=invblock[l];
		     if (!idxpos[m]) {
		        idx[cnt]=m;
			idxpos[m]=++cnt;
		     } /* end if */
		 } /* end for k */
		 j++;
	   } /* end while */
	   B.ia[i+1]=cnt;
	   nnz+=cnt;
	   for (k=0; k<cnt; k++) {
	       m=idx[k];
	       idxpos[m]=0;
	   } /* end for k */
	   cnt=0;  
       } /* end for i */

       /* change B.ia from nnz to pointer */
       B.ia[0]=0;
       for (i=0; i<mynblocks; i++) 
	   B.ia[i+1]+=B.ia[i];
       /* FORTRAN notation for companion matrix B for ILUPACK reordering routines */
       for (i=0; i<=mynblocks; i++) 
	   B.ia[i]++;
       
       B.ja=(integer *)malloc((size_t)nnz*sizeof(integer));
       B.a =(FLOAT *)  malloc((size_t)nnz*sizeof(FLOAT));
       B.nnz=nnz;
       /* second pass: detect pattern column by column */
       j=0;
       for (i=0; i<mynblocks; i++) {
	   /* scan block column i */
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         jj=p[j]; /* shortcut */
		 for (k=0; k<C->ncol[jj]; k++) {
		     ll=C->rowind[jj][k];
		     /* we will have q[l]=ll, so in total we are considering C(q[l],p[j]) */
		     l=invq[ll];
		     m=invblock[l];
		     if (!idxpos[m]) {
		        idx[cnt]=m;
			idxpos[m]=++cnt;
		     } /* end if */
		 } /* end for k */
		 j++;
	   } /* end while */
	   /* remember FORTRAN notation for B */
	   l=B.ia[i]-1;
	   for (k=0; k<cnt; k++) {
	       m=idx[k];
	       idxpos[m]=0;
	       B.ja[l] =m+1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	       B.a[l++]=1.0; /* dummy */
#else
	       B.a[l].r=1.0; /* dummy */
	       B.a[l].i=0.0; /* dummy */
	       l++;
#endif
	   } /* end for k */
	   cnt=0;  
       } /* end for i */  
       free(idx);
       free(idxpos);
#ifdef PRINT_INFO
       printf("SYMBILU: companion matrix for symmetric reorderings computed\n");fflush(stdout);
#endif

#ifdef PRINT_INFO2
       printf("SYMBILU: call (scaling+) reordering\n");fflush(stdout);
       for (i=0; i<B.nc; i++) {
	   printf("column %8ld\n",i+1);
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       k=B.ja[j-1];
	       printf("%8ld",k);
	   }
	   printf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       printf("%8.1le",B.a[j-1].r);
	   }
	   printf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       printf("%8.1le",B.a[j-1].i);
	   }
	   printf("\n");
       }
#endif

       /* select permutation and scaling driver by user request */
       nB=mynblocks;
       /* dummies for row scaling */
       pcolscale2=(FLOAT *)malloc((size_t)nB*sizeof(FLOAT));
       prowscale2=pcolscale2;
       /* block permutation */
       p2   =(integer *)malloc((size_t)nB*sizeof(integer));
       invq2=(integer *)malloc((size_t)nB*sizeof(integer));
       
#ifdef PRINT_INFO
       printf("SYMBILU: apply reordering to the compressed graph\n");fflush(stdout);
#endif
       // if (!strcmp("mtmetis",options.ordering)) {
       if (options.ordering==PERM_MTMETIS) {
	  ierr=MYSYMPERMMMTMETIS(B, prowscale2,pcolscale2, p2,invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: mtmetis performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       // else if (!strcmp("metisn",options.ordering)) {
       else if (options.ordering==PERM_METISN) {
	  ierr=MYSYMPERMMMETISN(B, prowscale2,pcolscale2, p2,invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: metisn performed\n");fflush(stdout);
#endif
       }      
       // else if (!strcmp("metise",options.ordering)) {
       else if (options.ordering==PERM_METISE) {
	  ierr=MYSYMPERMMMETISE(B, prowscale2,pcolscale2, p2,invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: metise performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       // else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
	  ierr=MYSYMPERMMAMD(B, prowscale2,pcolscale2, p2,invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: amd performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
	  ierr=MYSYMPERMMRCM(B, prowscale2,pcolscale2, p2,invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMBILU: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  for (i=0; i<nB; i++)
	      p2[i]=invq2[i]=i+1;
#ifdef PRINT_INFO
	  printf("SYMBILU: no symmetric reordering used\n");fflush(stdout);
#endif
       }
       free(pcolscale2);
      
    
#ifdef PRINT_INFO
       printf("SYMBILU: block (scaling +) reordering computed, ierr=%ld\n",ierr);fflush(stdout);
#endif
#ifdef PRINT_INFO0
       printf("block p as returned\n");
       for (i=0; i<B.nc; i++) {
	   printf("%8ld",p2[i]);
       }
       printf("\n");
       fflush(stdout);
       printf("block invq as returned \n");
       fflush(stdout);
       for (i=0; i<B.nc; i++) {
	   printf("%8ld",invq2[i]);
       }
       printf("\n");
       fflush(stdout);
#endif
       free(B.ia);
       free(B.ja);
       free(B.a);
#ifdef PRINT_CHECK
       for (i=0;i<nB; i++) {
	   if (p2[i]<=0 || p2[i]>nB) {
	      printf("SYMBILU ERROR, wrong permutation product, step %ld, p2[%ld]=%ld\n", i,i,p2[i]);
	      fflush(stdout);
	   }
	   if (invq2[i]<=0 || invq2[i]>nB) {
	      printf("SYMBILU ERROR, wrong permutation product, step %ld, invq2[%ld]=%ld\n", i,i,invq2[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
       for (i=0; i<nB; i++) {
	   if (check_buff[p2[i]-1]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
	      fflush(stdout);
	   }
	   check_buff[p2[i]-1]=-1;
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
       for (i=0; i<nB; i++) {
	   if (check_buff[invq2[i]-1]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, invq2[%ld]=%ld\n",i,i,invq2[i]);
	      fflush(stdout);
	   }
	   check_buff[invq2[i]-1]=-1;
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
#endif


       /* expand permutations from block to scalar */
       invq3=(integer *)malloc((size_t)n*sizeof(integer));
       p3   =(integer *)malloc((size_t)n*sizeof(integer));
       blockind=(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
       /* extract scalar starting positions from blocksize */
       blockind[0]=0;
       for (i=0; i<mynblocks; i++)
	   blockind[i+1]=blockind[i]+myblocksize[i];
       /* expand block p2 to scalar p3 */
       j=0;
       for (i=0; i<mynblocks; i++) {
	   k=p2[i]-1; /* FORTRAN to C notation */
	   for (l=blockind[k]; l<blockind[k+1]; l++)
	       p3[j++]=l;
       } /* end for i */

    
       /* blocksize is now incoroprated by the permutations */
       if (options.blocksize==NULL && !options.blocking_strategy) { 
	  /* no further a priori block requested, transfer block information */
	  /* use invblock as buffer */
	  for (i=0; i<mynblocks; i++) {
	      k=p2[i]-1; /* FORTRAN to C notation */
	      invblock[i]=myblocksize[k];
	  }
	  memcpy(myblocksize, invblock, mynblocks*sizeof(integer));
       }
#ifdef PRINT_CHECK
       if (mynblocks>0 && myblocksize!=NULL) {
	  k=0;
	  for (j=0; j<mynblocks; j++) {
	      k+=myblocksize[j];
	      if (myblocksize[j]<=0 || myblocksize[j]>n)
		 printf("SYMBILU ERROR, wrong block size, permuted cosine, step %ld, blocksize[%ld]=%ld\n", j,j,myblocksize[j]);
	      if (k>n)
		 printf("SYMBILU ERROR, cummulative block size exceeds system size, permuted cosine, step %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
	  }  
	  if (k!=n)
	     printf("SYMBILU ERROR, cummulative block does not match system size, permuted cosine, nblocks %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
#ifdef PRINT_CHECK0
	  printf("myblocksize[0,...,%ld]\n",mynblocks-1);fflush(stdout);
	  for (j=0; j<mynblocks; j++)
	      printf("%4ld",myblocksize[j]);
	  printf("\n");fflush(stdout);
#endif
       }
#endif

    
       /* compute block q2 from block invq2 and use p2 as buffer */
       for (i=0; i<mynblocks; i++)
	   p2[invq2[i]-1]=i;
       j=0;
       for (i=0; i<mynblocks; i++) {
	   k=p2[i]; /* now q2 (stored in p2) already uses C notation */
	   for (l=blockind[k]; l<blockind[k+1]; l++)
	       invblock[j++]=l; /* compute scalar q3 and use invblock as buffer */
       } /* end for i */
       /* compute scalar invq3 from scalar q3 (stored in invblock) */
       for (i=0; i<n; i++)
	   invq3[invblock[i]]=i;
#ifdef PRINT_CHECK
       for (i=0;i<n; i++) {
	   if (invblock[i]<0 || invblock[i]>=n) {
	      printf("SYMBILU ERROR, wrong permutation product, step %ld, invblock[%ld]=%ld\n", i,i,invblock[i]);
	      fflush(stdout);
	   }
	   if (invq3[i]<0 || invq3[i]>=n) {
	      printf("SYMBILU ERROR, wrong permutation product, step %ld, invq3[%ld]=%ld\n", i,i,invq3[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invblock[i]]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, invblock[%ld]=%ld\n",i,i,invblock[i]);
	      fflush(stdout);
	   }
	   check_buff[invblock[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq3[i]]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, invq3[%ld]=%ld\n",i,i,invq3[i]);
	      fflush(stdout);
	   }
	   check_buff[invq3[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif

#ifdef PRINT_INFO
       printf("SYMBILU: permutation from reordering of the companion matrix expanded to full size\n");fflush(stdout);
#endif
    
    
       /* block structures are not needed anymore */
       free(blockind);
       free(invq2);
       free(p2);
	  
       /* replace permutations by their product with the symmetric permutations */
       /* p <- p[p3] use invblock as buffer */
       for (i=0; i<n; i++)
	   invblock[i]=p[p3[i]];
       memcpy(p, invblock, n*sizeof(integer));
       /* invq <- invq3[invq] use invblock as buffer */
       for (i=0; i<n; i++)
	   invblock[i]=invq3[invq[i]];
       memcpy(invq, invblock, n*sizeof(integer));
#ifdef PRINT_CHECK0
       printf("p\n");
       for (j=0; j<n; j++)
	   printf("%6ld",p[j]);
       printf("\ninvq\n");
       fflush(stdout);
       for (j=0; j<n; j++)
	   printf("%6ld",invq[j]);
       printf("\n");
       fflush(stdout);
#endif
#ifdef PRINT_CHECK
       for (j=0; j<n; j++) {
	   if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n) {
	      printf("SYMBILU ERROR, wrong p/invq after reorering the compressed graph, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j],invq[p[j]]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[p[i]]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
	      fflush(stdout);
	   }
	   check_buff[p[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq[i]]) {
	      printf("SYMBILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
	      fflush(stdout);
	   }
	   check_buff[invq[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif

#ifdef PRINT_INFO
       printf("SYMBILU: permutation product built\n");fflush(stdout);
#endif
    
       free(invblock);
       free(p3);
       free(invq3);

    
#ifdef PRINT_INFO
       printf("SYMBILU: matrix for BILU prepared\n");fflush(stdout);
#endif
    } /* end if blocksize!=0 */
    /* pcolscale not required anymore since C is already properly scaled */
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
    free(colscale);
#endif
#ifdef _PROFILING_
    time_reordering+=omp_get_wtime()-timeBegin;
#endif
    

      


    
#ifdef PRINT_INFO
    printf("call SYMBILU, droptol=%8.1le\n", options.droptol); fflush(stdout);
#endif
   

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif   
    /* blocksize is now incorporated by the permutations */
    if (options.blocksize==NULL && !options.cosine) {
       myblocksize=FREE(myblocksize);
       mynblocks=0;
    }
    if (options.blocksize==NULL) { 
       if (options.blocking_strategy) {
	  myblocksize=FREE(myblocksize);
	  mynblocks=0;
       }
    }
    else {
       myblocksize=options.blocksize;
       mynblocks=options.nblocks;
    }
    // printf("ilu1t=%ld, mynblocks=%ld, blocksize=%ld\n",options.blocking_strategy, mynblocks,myblocksize); 
    ierr=MYSYMBILU(C, BL,BiD,BL, SL,SL, p,invq, determinant, isdefinite, myblocksize,mynblocks,
		   options.droptol,options.blocking_strategy,options.progressive_aggregation,options.perturbation,pivots,options.level_of_fill);
#ifdef PRINT_INFO
    printf("SYMBILU: bilu completed, ierr=%ld\n",ierr);fflush(stdout);
#endif
    if (options.blocksize==NULL && !options.blocking_strategy) {
       myblocksize=FREE(myblocksize);
       mynblocks=0;
    }
       
    if (ierr) {
       printf("SYMBILU: block ILU abnormally terminated");
       
       ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
       ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
       ilupackoptions.iaux =FREE(ilupackoptions.iaux);
       ilupackoptions.daux =FREE(ilupackoptions.daux);
       
       return (ierr);
    }
#ifdef _PROFILING_
    time_bilu=omp_get_wtime()-timeBegin;
#endif


#ifdef _PROFILING_
    time_total=omp_get_wtime()-time_total;
    printf("BILDL driver profiling summary\n");
    printf("  i) pre processing         %8.1le\n",time_pre_processing);
    printf(" ii) cosine-based blocking  %8.1le\n",time_cosine);
    printf("iii) reordering             %8.1le\n",time_reordering);
    printf(" iv) BILDL                  %8.1le\n",time_bilu);
    printf("Total BILDL driver time %8.1le\n\n",time_total);

    fflush(stdout);
#endif

#ifdef PRINT_CHECK
    free(check_buff);
#endif

    ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
    ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
    ilupackoptions.iaux =FREE(ilupackoptions.iaux);
    ilupackoptions.daux =FREE(ilupackoptions.daux);
       
    return (0);       

} // end bildldriver
