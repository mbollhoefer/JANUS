/* $Id: biludriver.c 7315 2021-05-28 21:00:20Z bolle $ 

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


// #define PRINT_INFO

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define USE_COSINE_BLOCKS 
#define TAU      0.8

#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

integer BILUDRIVER(SPARSEMATRIX *C,
		   SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT, integer *pivots,
		   REALS *SL, REALS *SR, integer *p, integer *invq,
                   doublecomplex *determinant,
		   JanusOptions options)
{  
    integer            i,j,k,l,ll,jj,m,n,*invq2,*p2,*invq3,*p3,r,s,mynblocks,*myblocksize=NULL,
                       nB=0, ierr,bi, *invblock, *blockind, nnz,
                       cnt, *idx, *idxpos;
    FLOAT              *pcolscale, *prowscale,
                       *pcolscale2, *prowscale2;
    size_t             mrows, ncols;
    ILUPACKPARAM       ilupackoptions;
    CSRMAT             A,B;
#ifdef _PROFILING_
  double timeBegin,
         time_total=0.0,
         time_bilu=0.0,
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
    
    /* copy input matrix to sparse row format */
    A.nc=A.nr=n;
    A.ia=(integer *)malloc((size_t)(n+1)*sizeof(integer));
    A.ja=(integer *)malloc((size_t)nnz     *sizeof(integer));
    A. a=(FLOAT *)  malloc((size_t)nnz     *sizeof(FLOAT));
    
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
    AMGINIT(&A,&ilupackoptions);
    if (options.symmetric_structure)
       ilupackoptions.ipar[6]|=SYMMETRIC_STRUCTURE;


    /* memory for initial scaling and permutation */
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    pcolscale=SL;
    prowscale=SR;
#else
    pcolscale=(FLOAT *)malloc((size_t)n*sizeof(FLOAT));
    prowscale=(FLOAT *)malloc((size_t)n*sizeof(FLOAT));
#endif
    nB=n;
#ifdef _PROFILING_
    time_pre_processing=omp_get_wtime()-timeBegin;
#endif

    
    /* scalar case, neither predefined blocks nor cosine strategy */
    if (!options.cosine && (options.nblocks==n||options.nblocks==0||options.blocksize==NULL)) {
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* remember that A is already prepared in FORTRAN-style indexing in order
	  to fit ILUPACK requirements */
       /* matching turned on */
       if (options.matching) {
#ifdef PRINT_INFO
	  printf("BILU: use maximum weight matching\n");fflush(stdout);
#endif
	  // if (!strcmp("mtmetis",options.ordering)) {
	  if (options.ordering==PERM_MTMETIS) {
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64MTMETIS(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGMTMETIS(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + MT-METIS performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  // else if (!strcmp("metisn",options.ordering)) {
	  else if (options.ordering==PERM_METISN) {
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64METISN(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGMETISN(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + METIS-N performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("metise",options.ordering)) {
	  else if (options.ordering==PERM_METISE) {
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64METISE(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGMETISE(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + METIS-E performed\n");fflush(stdout);
#endif
	  }
#endif // _USE_METIS4_
	  // else if (!strcmp("amd",options.ordering)) {
	  else if (options.ordering==PERM_AMD) {
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64AMD(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGAMD(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + AMD performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("rcm",options.ordering)) {
	  else if (options.ordering==PERM_RCM) {
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64RCM(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGRCM(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + RCM performed\n");fflush(stdout);
#endif
	  }
	  else { /* none */
#ifdef _MC64_MATCHING_    
	     ierr=PERMMC64NULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	     ierr=PERMMATCHINGNULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: just MC64 + scaling performed\n");fflush(stdout);
#endif
	  }
       }
       else { /* no matching */
	  
	  // if (!strcmp("mtmetis",options.ordering)) {
	  if (options.ordering==PERM_MTMETIS) {
	     ierr=PERMMTMETIS(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: MT-METIS performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  // else if (!strcmp("metisn",options.ordering)) {
	  else if (options.ordering==PERM_METISN) {
	     ierr=PERMMETISN(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: METIS-N performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("metise",options.ordering)) {
	  else if (options.ordering==PERM_METISE) {
	     ierr=PERMMETISE(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: METIS-E performed\n");fflush(stdout);
#endif
	  }
#endif //_USE_METIS4_
	  // else if (!strcmp("amd",options.ordering)) {
	  else if (options.ordering==PERM_AMD) {
	     ierr=PERMAMD(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: AMD performed\n");fflush(stdout);
#endif
	  }
	  // else if (!strcmp("rcm",options.ordering)) {
	  else if (options.ordering==PERM_RCM) {
	     ierr=PERMRCM(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: RCM performed\n");fflush(stdout);
#endif
	  }
	  else { /* none */
	     ierr=PERMNULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	     printf("BILU: just scaling performed\n");fflush(stdout);
#endif
	  }
       }
       /* switch permutations to C-style */
       for (i=0; i<n; i++) {
	   p[i]--;
	   invq[i]--;
       } /* end for i */
#ifdef _PROFILING_
       time_reordering=omp_get_wtime()-timeBegin;
#endif

#ifdef PRINT_INFO
       printf("BILU: C set up\n");fflush(stdout);
#endif

      

#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
       for (i=0; i<n; i++) {
	   SR[i]=prowscale[i].r;
	   SL[i]=pcolscale[i].r;
       } // end for i
       free(prowscale);
       free(pcolscale);
#endif
       ierr=BILU(C, BL,BiD,BUT, SL,SR, p,invq, determinant, NULL,0,
		 options.droptol,options.blocking_strategy,options.progressive_aggregation,options.perturbation,options.symmetric_structure,pivots,options.level_of_fill);
       
       /* matrix A is not longer required */
       free(A.ia);  
       free(A.ja);
       free(A.a);
       
       if (ierr) {
	  printf("BILU: block ILU abnormally terminated\n");
	  
	  ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	  ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	  ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	  ilupackoptions.daux =FREE(ilupackoptions.daux);
	  
	  return (ierr);
       }
#ifdef PRINT_INFO
       else {
	  printf("BILU: block ILU computation successfully completed\n");fflush(stdout);
       }
#endif

#ifdef _PROFILING_
       time_bilu=omp_get_wtime()-timeBegin;
#endif
    } /* end if neither cosine nor initial blocks */
    else {
       /* block case, predefined blocks or cosine strategy */

#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* compute maximum weight matching as well as scalings */
       if (options.matching) {
#ifdef PRINT_INFO
	  printf("BILU: use maximum weight matching\n");fflush(stdout);
#endif
#ifdef _MC64_MATCHING_    
	  ierr=PERMMC64NULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#else
	  ierr=PERMMATCHINGNULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
#endif
#ifdef PRINT_INFO
	     printf("BILU: MC64 + scaling performed\n");fflush(stdout);
#endif
       }
       else { /* extract the diagonal scaling only */
#ifdef PRINT_INFO
	  printf("BILU: no maximum weight matching, just scaling\n");fflush(stdout);
#endif
	  ierr=PERMNULL(A, prowscale,pcolscale, p,invq, &nB, &ilupackoptions);
       }
#ifdef PRINT_INFO
       printf("BILU: scaling and permutation prior to block partioning computed, ierr=%d\n",ierr);fflush(stdout);
#endif
       /* switch permutations to C-style */
       for (i=0; i<n; i++) {
	   p[i]--;
	   invq[i]--;
       } /* end for i */

       /* switch matrix A back to C notation and adapt data structures for bilu (and cosine_blocks) */
       for (i=0; i<=n; i++)
	   A.ia[i]--;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]--;
#ifdef _PROFILING_
       time_pre_processing+=omp_get_wtime()-timeBegin;
#endif

    
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* if the user does not predefine an initial block structure we do Saad's cosine blocking */
       if (options.blocksize==NULL) {
	  myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
	  mynblocks=0;
	  /* previously computed permutations p and invq will be respected and updated */
#ifdef USE_COSINE_BLOCKS       
	  if (options.cosine) {
#ifdef PRINT_INFO
	     printf("BILU: call modified cosine-based blocking\n");fflush(stdout);
#endif
	     MCOSINE_BLOCKS(C, p,invq, myblocksize, &mynblocks, TAU);
#ifdef PRINT_INFO
	     printf("BILU: modified cosine-based blocking applied\n");fflush(stdout);
#endif
	  }
	  else {
	     for (i=0; i<n; i++)
	         myblocksize[i]=1;
	     mynblocks=n;
	  } /* end if-else */
#else
	  for (i=0; i<n; i++)
	      myblocksize[i]=1;
	  mynblocks=n;
#endif
       }
       else {
	  /* remove possible non-symmetric reordering from matchings */
	  for (i=0; i<n; i++)
	      p[i]=invq[i]=i;
	  myblocksize=options.blocksize;
	  mynblocks=options.nblocks;
#ifdef PRINT_INFO
	  printf("BILU: initial permutation turned off\n");fflush(stdout);
#endif
       }
#ifdef _PROFILING_
       time_cosine=omp_get_wtime()-timeBegin;
#endif

       
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* compute a companion matrix B of A(q,p) compressing the associated blocks */
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
       B.ia=(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
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
		 for (k=A.ia[jj]; k<A.ia[jj+1]; k++) {
		     ll=A.ja[k];
		     /* we will have q[l]=ll, so in total we are considering A(q[l],p[j]) */
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
		 for (k=A.ia[jj]; k<A.ia[jj+1]; k++) {
		     ll=A.ja[k];
		     /* we will have q[l]=ll, so in total we are considering A(q[l],p[j]) */
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
	       B.a[l].r  =1.0; /* dummy */
	       B.a[l++].i=0.0; /* dummy */
#endif
	   } /* end for k */
	   cnt=0;  
       } /* end for i */  
       free(idx);
       free(idxpos);
#ifdef PRINT_INFO
       printf("BILU: companion matrix computed\n");fflush(stdout);
#endif


       /* select permutation and scaling driver by user request */
       nB=mynblocks;
       /* block permutation */
       p2   =(integer *)malloc((size_t)nB*sizeof(integer));
       invq2=(integer *)malloc((size_t)nB*sizeof(integer));
       /* dummies for row scaling */
       pcolscale2=(FLOAT *)malloc((size_t)nB*sizeof(FLOAT));
       prowscale2=(FLOAT *)malloc((size_t)nB*sizeof(FLOAT));
       // if (!strcmp("mtmetis",options.ordering)) {
       if (options.ordering==PERM_MTMETIS) {
	  ierr=PERMMTMETIS(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("BILU: MT-METIS performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       // else if (!strcmp("metisn",options.ordering)) {
       else if (options.ordering==PERM_METISN) {
	  ierr=PERMMETISN(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("BILU: METIS-N performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("metise",options.ordering)) {
       else if (options.ordering==PERM_METISE) {
	  ierr=PERMMETISE(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("BILU: METIS-E performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       //else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
	  ierr=PERMAMD(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("BILU: AMD performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
	  ierr=PERMRCM(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("BILU: RCM performed\n");fflush(stdout);
#endif
       }
       else { /* none */
          for (i=0; i<nB; i++)
	      p2[i]=invq2[i]=i+1;
#ifdef PRINT_INFO
	  printf("BILU: no symmetric reordering used\n");fflush(stdout);
#endif
       }
       free(prowscale2);
       free(pcolscale2);
    
#ifdef PRINT_INFO
       printf("BILU: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
       free(B.ia);
       free(B.ja);
       free(B.a);


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

       /* block structures are not needed anymore */
       free(blockind);
       free(invq2);
       free(p2);
#ifdef PRINT_INFO
       printf("BILU: permutation expanded\n");fflush(stdout);
#endif
    
       /* replace permutations by their product with the symmetric permutations */
       /* p <- p[p3] use invblock as buffer */
       for (i=0; i<n; i++)
           invblock[i]=p[p3[i]];
       memcpy(p, invblock, n*sizeof(integer));
       /* invq <- invq3[invq] use invblock as buffer */
       for (i=0; i<n; i++)
	   invblock[i]=invq3[invq[i]];
       memcpy(invq, invblock, n*sizeof(integer));
    
       free(invblock);
       free(p3);
       free(invq3);
#ifdef PRINT_INFO
       printf("BILU: permutation product built\n");fflush(stdout);
#endif

    

#ifdef _PROFILING_
       time_reordering=omp_get_wtime()-timeBegin;
#endif
    


    
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif   
#ifdef PRINT_INFO
       printf("call BILU, droptol=%8.1le\n", options.droptol); fflush(stdout);
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
                
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
       for (i=0; i<n; i++) {
	   SR[i]=prowscale[i].r;
	   SL[i]=pcolscale[i].r;
       } // end for i
       free(prowscale);
       free(pcolscale);
#endif
       ierr=BILU(C, BL,BiD,BUT, SL,SR, p,invq, determinant, myblocksize,mynblocks,
		 options.droptol,options.blocking_strategy,options.progressive_aggregation,options.perturbation,options.symmetric_structure,pivots,options.level_of_fill);
#ifdef PRINT_INFO
       printf("BILU: bilu completed, ierr=%d\n",ierr);fflush(stdout);
#endif
       if (options.blocksize==NULL && !options.blocking_strategy) {
	  myblocksize=FREE(myblocksize);
	  mynblocks=0;
       }
       /* matrix A is not longer required */
       free(A.ia);  
       free(A.ja);
       free(A.a);

       if (ierr) {
	  printf("BILUdriver: block ILU abnormally terminated\n");
	  
	  ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	  ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	  ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	  ilupackoptions.daux =FREE(ilupackoptions.daux);
	  
	  return (ierr);
       }
#ifdef PRINT_INFO
       else {
	  printf("BILU: block ILU computation successfully completed\n");fflush(stdout);
       }
#endif
#ifdef _PROFILING_
       time_bilu=omp_get_wtime()-timeBegin;
#endif
    } /* end-if-else scalar vs. block case */

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif   

    


#ifdef _PROFILING_
    time_total=omp_get_wtime()-time_total;
    printf("BILU driver profiling summary\n");
    printf("  i) pre processing         %8.1le\n",time_pre_processing);
    printf(" ii) cosine-based blocking  %8.1le\n",time_cosine);
    printf("iii) reordering             %8.1le\n",time_reordering);
    printf(" iv) BILU                   %8.1le\n",time_bilu);
    printf("Total BILU driver time %8.1le\n\n",time_total);

    fflush(stdout);
#endif


    ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
    ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
    ilupackoptions.iaux =FREE(ilupackoptions.iaux);
    ilupackoptions.daux =FREE(ilupackoptions.daux);
    
    return (0);       

}
