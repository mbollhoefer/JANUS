/* $Id: bicholdriver.c 7315 2021-05-28 21:00:20Z bolle $ 

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

#define USE_COSINE_BLOCKS 
#define TAU      0.8
// #define PRINT_INFO

#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

integer SPDBILUDRIVER(SPARSEMATRIX *C,
		      SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT, 
		      REALS *SL, REALS *SR, integer *p, integer *invq, 
		      doublecomplex *determinant,
		      JanusOptions options)
{  
    integer            i,j,k,l,ll,jj,m,n,*invq2,*p2,*invq3,*p3,r,s,mynblocks,
                       *myblocksize=NULL, nB=0, ierr,bi, *invblock, *blockind, nnz,
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
    
    /* so far no reordering */
    for (i=0; i<n; i++)
        p[i]=invq[i]=i;
    
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
    SPDAMGINIT(&A,&ilupackoptions);
#ifdef _PROFILING_
    time_pre_processing=omp_get_wtime()-timeBegin;
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
#ifdef PRINT_INFO
       printf("bicholdriver: call modified cosine-based blocking\n");fflush(stdout);
#endif
       if (options.cosine) {
#ifdef PRINT_INFO
	  printf("Call MCosine_SBlocks\n");
#endif
	  MCOSINE_SBLOCKS(C, p,invq, myblocksize, &mynblocks, TAU);
       }
       else {
	  for (i=0; i<n; i++)
	      myblocksize[i]=1;
	  mynblocks=n;
       } /* end if-else */
#ifdef PRINT_INFO
       printf("bicholdriver: modified cosine-based blocking applied\n");fflush(stdout);
#endif
#else
       for (i=0; i<n; i++)
	   myblocksize[i]=1;
       mynblocks=n;
#endif
    }
    else {
       myblocksize=options.blocksize;
       mynblocks=options.nblocks;
#ifdef PRINT_INFO
       printf("bicholdriver: initial permutation turned off\n");fflush(stdout);
#endif
    }
#ifdef _PROFILING_
    time_cosine=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
    /* row scaling vector */
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    pcolscale=prowscale=SL;
#else
    pcolscale=(FLOAT *)malloc((size_t)n*sizeof(FLOAT));
    prowscale=pcolscale;
#endif

    /* cosine not applied or cosine did not find dense blocks */
    if (mynblocks==n) {
       /* select permutation and scaling driver by user request */
       // if (!strcmp("mtmetis",options.ordering)) {
       if (options.ordering==PERM_MTMETIS) {
#ifdef PRINT_INFO
	  printf("bicholdriver: perform mtmetis\n");fflush(stdout);
#endif
          ierr=SPDPERMMTMETIS(A, prowscale, pcolscale, p, invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: mtmetis performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       // else if (!strcmp("metisn",options.ordering)) {
       else if (options.ordering==PERM_METISN) {
	  ierr=SPDPERMMETISN(A, prowscale, pcolscale, p, invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: metisn performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("metise",options.ordering)) {
       else if (options.ordering==PERM_METISE) {
          ierr=SPDPERMMETISE(A, prowscale, pcolscale, p, invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: metise performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       // else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
          ierr=SPDPERMAMD(A, prowscale, pcolscale, p, invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: amd performed\n");fflush(stdout);
#endif
       }
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
	  ierr=SPDPERMRCM(A, prowscale, pcolscale, p, invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  ierr=SPDPERMNULL(A, prowscale,pcolscale, p,invq, &n, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: no maximum weight matching, just scaling\n");fflush(stdout);
#endif
       }

       if (ierr) {
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
       
       /* switch back to C to style */
       for (i=0; i<=n; i++)
	   A.ia[i]--;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]--;
       for (i=0; i<n; i++) {
	   p[i]--;
	   invq[i]--;
       } /* end for i */
       

#ifdef PRINT_INFO
       printf("bicholdriver: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
    
#ifdef PRINT_INFO
       printf("bicholdriver: permutation product built\n");fflush(stdout);
#endif

    
#ifdef PRINT_INFO
       printf("bicholdriver: matrix for BICHOL prepared\n");fflush(stdout);
#endif
    } // end if mynblocks=n
    else { /* cosine successful or user-defined blocking */
       /* permutation, right now ignored, but later used for companion matrix */
       p2   =(integer *)malloc((size_t)n*sizeof(integer));
       invq2=(integer *)malloc((size_t)n*sizeof(integer));
       ierr=SPDPERMNULL(A, prowscale,pcolscale, p2,invq2, &n, &ilupackoptions);

       if (ierr) {
	  myblocksize=FREE(myblocksize);
	  free(p2);
	  free(invq2);
	  free(A.ia);
	  free(A.ja);
	  free(A.a);
	  ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	  ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	  ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	  ilupackoptions.daux =FREE(ilupackoptions.daux);
	  return ierr;
       } // end if ierr
       
       /* switch back to C to style */
       for (i=0; i<=n; i++)
	   A.ia[i]--;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]--;
       
       /* compute a companion matrix B of A(p,p) compressing the associated blocks */
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
       printf("bicholdriver: companion matrix computed\n");fflush(stdout);
#endif


       /* select permutation and scaling driver by user request */
       nB=mynblocks;
       /* dummies for row scaling */
       pcolscale2=(FLOAT *)malloc((size_t)nB*sizeof(FLOAT));
       prowscale2=pcolscale2;
       // if (!strcmp("mtmetis",options.ordering)) {
       if (options.ordering==PERM_MTMETIS) {
          ierr=SPDPERMMTMETIS(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: mtmetis performed\n");fflush(stdout);
#endif
       } // end if-else if mtmetis
#ifdef _USE_METIS4_
       // else if (!strcmp("metisn",options.ordering)) {
       else if (options.ordering==PERM_METISN) {
	  ierr=SPDPERMMETISN(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: metisn performed\n");fflush(stdout);
#endif
       } // end if metisn
       // else if (!strcmp("metise",options.ordering)) {
       else if (options.ordering==PERM_METISE) {
          ierr=SPDPERMMETISE(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: metise performed\n");fflush(stdout);
#endif
       } // end if-else if metise
#endif // _USE_METIS4_
       // else if (!strcmp("amd",options.ordering)) {
       else if (options.ordering==PERM_AMD) {
          ierr=SPDPERMAMD(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: amd performed\n");fflush(stdout);
#endif
       } // end if-else if amd
       // else if (!strcmp("rcm",options.ordering)) {
       else if (options.ordering==PERM_RCM) {
	  ierr=SPDPERMRCM(B, prowscale2, pcolscale2, p2, invq2, &nB, &ilupackoptions);
#ifdef PRINT_INFO
	  printf("bicholdriver: rcm performed\n");fflush(stdout);
#endif
       } // end if-else if rcm
       else { /* none */
	  for (i=0; i<nB; i++)
	      p2[i]=invq2[i]=i+1;
#ifdef PRINT_INFO
	  printf("bicholdriver: no symmetric reordering used\n");fflush(stdout);
#endif
       } // end if-else if-else mtmetis
       free(pcolscale2);
    
#ifdef PRINT_INFO
       printf("bicholdriver: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
       free(B.ia);
       free(B.ja);
       free(B.a);

       if (ierr) {
	  
	  myblocksize=FREE(myblocksize);
	  free(invblock);
	  free(p2);
	  free(invq2);
	  free(A.ia);
	  free(A.ja);
	  free(A.a);
	  ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
	  ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
	  ilupackoptions.iaux =FREE(ilupackoptions.iaux);
	  ilupackoptions.daux =FREE(ilupackoptions.daux);
	  return ierr;
       } // end if ierr


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
	  } // end for i
	  memcpy(myblocksize, invblock, mynblocks*sizeof(integer));
       } // end options.blocksize=0

    
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

#ifdef PRINT_INFO
       printf("bicholdriver: permutation expanded\n");fflush(stdout);
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
    
       free(invblock);
       free(p3);
       free(invq3);
#ifdef PRINT_INFO
       printf("bicholdriver: permutation product built\n");fflush(stdout);
#endif

    
#ifdef PRINT_INFO
       printf("bicholdriver: matrix for BICHOL prepared\n");fflush(stdout);
#endif
    } // if-else mynblocks=n
#ifdef _PROFILING_
    time_reordering=omp_get_wtime()-timeBegin;
#endif
    

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
    for (i=0; i<n; i++) {
        SL[i]=pcolscale[i].r;
    } // end for i
    free(pcolscale);
#endif

#ifdef PRINT_INFO
    printf("call BICHOL, droptol=%8.1le\n", options.droptol); fflush(stdout);
#endif


    if (options.blocksize==NULL && !options.cosine) {
       myblocksize=FREE(myblocksize);
       mynblocks=0;
    } // end if blocksize=0					   
    if (options.blocksize==NULL) { 
       if (options.blocking_strategy) {
	  myblocksize=FREE(myblocksize);
	  mynblocks=0;
       } // end if blocking_strategy
    } // end if blocksize=0					   
    else {
       myblocksize=options.blocksize;
       mynblocks=options.nblocks;
    } // end if-else blocksize=0					   
    /* printf("blocking_strategy=%ld, mynblocks=%ld, blocksize=%ld\n",options.blocking_strategy, mynblocks,myblocksize);fflush(stdout); */
    ierr=BILDLSPD(C, BL,BiD,BL, SL,SL, p,invq, determinant, myblocksize,mynblocks,
		  options.droptol,options.blocking_strategy,options.progressive_aggregation,options.perturbation,options.invert_blocks,options.level_of_fill);
#ifdef PRINT_INFO
    printf("bicholdriver: bilu completed, ierr=%d\n",ierr);fflush(stdout);
#endif
    myblocksize=FREE(myblocksize);
    mynblocks=0;
    /* matrix A is not longer required */
    free(A.ia);  
    free(A.ja);
    free(A.a);

    if (ierr) {
       printf("bicholdriver: block ILU abnormally terminated\n") ;

       ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
       ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
       ilupackoptions.iaux =FREE(ilupackoptions.iaux);
       ilupackoptions.daux =FREE(ilupackoptions.daux);
       
       return (ierr);
    } // end if ierr
#ifdef _PROFILING_
    time_bilu=omp_get_wtime()-timeBegin;
#endif


#ifdef _PROFILING_
    time_total=omp_get_wtime()-time_total;
    printf("BICHOL driver profiling summary\n");
    printf("  i) pre processing         %8.1le\n",time_pre_processing);
    printf(" ii) cosine-based blocking  %8.1le\n",time_cosine);
    printf("iii) reordering             %8.1le\n",time_reordering);
    printf(" iv) BICHOL                 %8.1le\n",time_bilu);
    printf("Total BICHOL driver time %8.1le\n\n",time_total);

    fflush(stdout);
#endif



    ilupackoptions.ibuff=FREE(ilupackoptions.ibuff);
    ilupackoptions.dbuff=FREE(ilupackoptions.dbuff);
    ilupackoptions.iaux =FREE(ilupackoptions.iaux);
    ilupackoptions.daux =FREE(ilupackoptions.daux);
       
    return (0);       

} // end bicholdriver
