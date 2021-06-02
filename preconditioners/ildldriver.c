/* $Id: ildldriver.c 4004 2018-02-14 18:51:05Z bolle $ */
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
#define MYSYMILU       SSYMilu
#define MYILDLDRIVER   SSYMilu_driver
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
#define MYSYMILU       DSYMilu
#define MYILDLDRIVER   DSYMilu_driver
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
#define MYSYMILU       CSYMilu
#define MYILDLDRIVER   CSYMilu_driver
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
#define MYSYMILU       ZSYMilu
#define MYILDLDRIVER   ZSYMilu_driver
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
#define MYSYMILU       SSYMilu
#define MYILDLDRIVER   SSYMilu_driver
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
#define MYSYMILU       DSYMilu
#define MYILDLDRIVER   DSYMilu_driver
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
#define MYSYMILU       CHERilu
#define MYILDLDRIVER   CHERilu_driver
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
#define MYSYMILU       ZHERilu
#define MYILDLDRIVER   ZHERilu_driver
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

integer MYILDLDRIVER(SPARSEMATRIX *C,
		     SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT,
		     REALS *SL, REALS *SR, integer *p, integer *invq, 
		     doublecomplex *determinant,
		     integer *isdefinite,
		     JanusOptions options)
{  
    integer            i,j,k,l,ii,ll,jj,m,n,r,s,
                       nB=0, ierr,bi, nnz, 
                       cnt, *idx, *idxpos;
    FLOAT              *colscale, *rowscale;
    size_t             mrows, ncols;
    ILUPACKPARAM       ilupackoptions;
    CSRMAT             A;
#ifdef PRINT_CHECK
    integer *check_buff;
#endif
#ifdef _PROFILING_
    double timeBegin,
           time_total=0.0,
           time_ilu=0.0,
           time_pre_processing=0.0,
           time_cosine=0.0,
           time_reordering=0.0;

    timeBegin=time_total=omp_get_wtime();
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

    /* row scaling vector */
    colscale=(FLOAT *)malloc((size_t)n*sizeof(FLOAT));
    rowscale=colscale;
    nB=n;
#ifdef _PROFILING_
    time_pre_processing+=omp_get_wtime()-timeBegin;
#endif

    
#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
    if (options.matching) { /* only matching */
       // if (!strcmp("metise",options.ordering)) {
       if (options.ordering==PERM_METISE) {
#ifdef _MC64_MATCHING_
	  ierr=MYSYMPERMMC64METISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	  ierr=MYSYMPERMMATCHINGMETISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	  printf("SYMILU: matching+metise performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("mtmetis",options.ordering)) {
       else if (options.ordering==PERM_MTMETIS) {
#ifdef _MC64_MATCHING_
	  ierr=MYSYMPERMMC64MTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	  ierr=MYSYMPERMMATCHINGMTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	  printf("SYMILU: matching+mtmetis performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
#ifdef _MC64_MATCHING_
	  ierr=MYSYMPERMMC64AMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	  ierr=MYSYMPERMMATCHINGAMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	  printf("SYMILU: matching+amd performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
#ifdef _MC64_MATCHING_
	  ierr=MYSYMPERMMC64RCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	  ierr=MYSYMPERMMATCHINGRCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	  printf("SYMILU: matching+rcm performed\n");fflush(stdout);
#endif
       }
       else { /* METIS nested dissection by nodes */
#ifdef _MC64_MATCHING_
	  ierr=MYSYMPERMMC64METISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#else
	  ierr=MYSYMPERMMATCHINGMETISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#endif	     
#ifdef PRINT_INFO
	  printf("SYMILU: matching+metisn performed\n");fflush(stdout);
#endif
       }
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
#ifdef PRINT_CHECK
       for (i=0;i<n; i++) {
	   if (p[i]<0 || p[i]>=n) {
	      printf("SYMILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
	      fflush(stdout);
	   }
	   if (invq[i]<0 || invq[i]>=n) {
	      printf("SYMILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[p[i]]) {
	      printf("SYMILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
	      fflush(stdout);
	   }
	   check_buff[p[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq[i]]) {
	      printf("SYMILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
	      fflush(stdout);
	   }
	   check_buff[invq[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif

       
       /* we are ready to apply ILDL */
    } /* end if only matching */
    else { /* no matching */
       /* do scaling and simple reordering without asking for matchings */
       // if (!strcmp("metisn",options.ordering)) {
       if (options.ordering==PERM_METISN) {
	  ierr=MYSYMPERMMMETISN(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: metisn performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("metise",options.ordering)) {
       else if (options.ordering==PERM_METISE) {
	  ierr=MYSYMPERMMMETISE(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: metise performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("mtmetis",options.ordering)) {
       else if (options.ordering==PERM_MTMETIS) {
	  ierr=MYSYMPERMMMTMETIS(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: mtmetis performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
	  ierr=MYSYMPERMMAMD(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: amd performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
	  ierr=MYSYMPERMMRCM(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  ierr=MYSYMPERMMNULL(A, rowscale,colscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("SYMILU: no symmetric reordering used\n");fflush(stdout);
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
	      printf("SYMILU ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
	      fflush(stdout);
	   }
	   if (invq[i]<0 || invq[i]>=n) {
	      printf("SYMILU ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[p[i]]) {
	      printf("SYMILU ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
	      fflush(stdout);
	   }
	   check_buff[p[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq[i]]) {
	      printf("SYMILU ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
	      fflush(stdout);
	   }
	   check_buff[invq[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif
       /* we are ready to apply ILDL */
    }

    /* matrix A is not longer required */
    free(colscale);  
    free(A.ia);  
    free(A.ja);
    free(A.a);


#ifdef _PROFILING_
    time_reordering+=omp_get_wtime()-timeBegin;
#endif
    
      

    
#ifdef PRINT_INFO
    printf("call SYMILU, droptol=%8.1le\n", options.droptol); fflush(stdout);
#endif
   

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif   
    ierr=MYSYMILU(C, BL,BiD,BL, SL,SL, p,invq, determinant, isdefinite, options.droptol);
#ifdef PRINT_INFO
    printf("SYMILU: ilu completed, ierr=%ld\n",ierr);fflush(stdout);
#endif
       
    if (ierr) {
       printf("SYMILU: block ILU abnormally terminated");
       return (ierr);
    }
#ifdef _PROFILING_
    time_ilu=omp_get_wtime()-timeBegin;
#endif


#ifdef _PROFILING_
    time_total=omp_get_wtime()-time_total;
    printf("ILDL driver profiling summary\n");
    printf("  i) pre processing         %8.1le\n",time_pre_processing);
    printf(" ii) reordering             %8.1le\n",time_reordering);
    printf("iii) ILDL                   %8.1le\n",time_ilu);
    printf("Total ILDL driver time %8.1le\n\n",time_total);

    fflush(stdout);
#endif

#ifdef PRINT_CHECK
    free(check_buff);
#endif

    return (0);       

} // end ildldriver
