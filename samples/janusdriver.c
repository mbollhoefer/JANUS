#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"

#include <janus.h>


#define PERTURB    1L
#define RES_TOL    1e-6
#define MAX_IT     5000L
#define K_RESTART    30L

#define MAX(A,B)        (((A)>(B))?(A):(B))


#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif


int main(int argc, char **argv)
{

    // matrix market variables to read in a sparse matrix
    MM_INT ret_code;
    MM_typecode matcode;
    FILE *f;
    MM_INT M, N, nz;   
    MM_INT *I, *J;
    double *valr,*vali;

    // major Janus variables
    SparseMatrix    A;
    JanusOptions options;
    JanusPrec PREC;
    void *sol, *rhs;

    // some minor variables
    long int nblocks;
    long int i,j,k,l,m,n,ierr,mx, iter;
    doublecomplex *pz, **pzv;
    double        *px, **pdv, v,
                  fill_bilu,
                  mu, stddev;
    long int CS, PA, LEVEL_OF_FILL;
    double DROP_TOL;

#ifdef _PROFILING_
    double timeBegin,
           time_total=0.0,
           time_bilu=0.0,
           time_krylov=0.0;
#endif


    // read in matrix in matrix market format
    if (argc<6) {
       fprintf(stderr, "Usage: %s [matrix-market-filename] [COSINE 1/0] [LEVEL_OF_FILL] [PROGR. AGGR. 1/0] [DROP_TOL]\n", argv[0]);
       exit(1);
    }
    else { 
       if ((f=fopen(argv[1], "r"))==NULL) 
	  exit(1);
    }
    if (mm_read_banner(f, &matcode)!=0){
       printf("Could not process Matrix Market banner.\n");
       exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode)) {
       printf("Sorry, this application does not support ");
       printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
       exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code=mm_read_mtx_crd_size(f, &M,&N,&nz))!=0)
       exit(1);
    if (M!=N) {
       printf("matrix must be square\n");
       exit(1);
    }

    /* reserve memory for matrices */
    I  =(MM_INT *)malloc((size_t)nz*sizeof(MM_INT));
    J  =(MM_INT *)malloc((size_t)nz*sizeof(MM_INT));
    if (mm_is_real(matcode)) {
       valr=(double *)malloc((size_t)nz*sizeof(double));
    }
    else {
       valr=(double *)malloc((size_t)nz*sizeof(double));
       vali=(double *)malloc((size_t)nz*sizeof(double));
    }

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++) {
        if (mm_is_real(matcode)) 
           fscanf(f, "%ld %ld %lg\n", &I[i], &J[i], &valr[i]);
	else
           fscanf(f, "%ld %ld %lg %lg\n", &I[i], &J[i], &valr[i],&vali[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f!=stdin) fclose(f);



    // convert matrix to my own sparse matrix format
    SparseInit(&A);
    n=N;
    A.nr=M; A.nc=N;
    A.nnz=nz;

    if (mm_is_real(matcode)) {
       A.isreal=1;
       A.val   =(void **) malloc((size_t)n*sizeof(double *));
       pdv=(double **) A.val;
    }
    else {
       A.isreal=0;
       A.val   =(void **) malloc((size_t)n*sizeof(doublecomplex *));
       pzv=(doublecomplex **)A.val;
    }

    // !!!!! ------------------------------------------ !!!!!
    // !!!!! ------------------------------------------ !!!!!
    // !!!!! ------------------------------------------ !!!!!
    // change this property to 1 if you know that your matrix
    // is (Hermitian) AND POSITIVE DEFINITE, otherwise use 0!
    A.isdefinite=0;
    // !!!!! ------------------------------------------ !!!!!
    // !!!!! ------------------------------------------ !!!!!
    // !!!!! ------------------------------------------ !!!!!

    if (mm_is_symmetric(matcode)) 
       A.issymmetric=1;
    else
       A.issymmetric=0;
    if (mm_is_hermitian(matcode)) 
       A.ishermitian=1;
    else
       A.ishermitian=0;
    
    A.ncol  =(integer *) calloc((size_t)n,sizeof(integer));
    A.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
    
    // find out number of entries per column
    for (k=0; k<nz; k++) {
        j=J[k];
	A.ncol[j]++;
    } // end for k
    // allocate memory for each column
    for (j=0; j<n; j++) {
        k=A.ncol[j];
	A.rowind[j]=(integer *)malloc((size_t)k*sizeof(integer));
	if (A.isreal)
	   pdv[j]=(double *) malloc((size_t)k*sizeof(double));
	else
	   pzv[j]=(doublecomplex *) malloc((size_t)k*sizeof(doublecomplex));
	A.ncol[j]=0;
    } // end for k

    // copy matrix
    for (k=0; k<nz; k++) {
        i=I[k];
        j=J[k];
	l=A.ncol[j];
	A.rowind[j][l]=i;
	A.ncol[j]=l+1;
	if (A.isreal)
	   pdv[j][l]=valr[k];
	else {
	   pzv[j][l].r=valr[k];
	   pzv[j][l].i=vali[k];
	}
    } // end for k
    free(I);
    free(J);
    if (A.isreal)
       free(valr);
    else {
       free(valr);
       free(vali);
    }
    if (A.ishermitian||A.issymmetric)
       if (A.isdefinite)
	  printf("apply BICHOL to %s, n=%ld,nz=%ld, rel. fill=%8.1le\n",argv[1],n,nz,(double)nz/n);
       else
	  printf("apply BILDL  to %s, n=%ld,nz=%ld, rel. fill=%8.1le\n",argv[1],n,nz,(double)nz/n);
    else
       printf("apply BILU   to %s, n=%ld,nz=%ld, rel. fill=%8.1le\n",argv[1],n,nz,(double)nz/n);


    // set up for Janus once we have read the matrix
#ifdef _PROFILING_
    time_total=omp_get_wtime();
#endif
    // allocate memory for right hand side and solution
    if (A.isreal) {
       sol=(void *) malloc((size_t)n*sizeof(double));
       rhs=(void *) malloc((size_t)n*sizeof(double));
       px =(double *) sol;
    }
    else {
       sol=(void *) malloc((size_t)n*sizeof(doublecomplex));
       rhs=(void *) malloc((size_t)n*sizeof(doublecomplex));
       pz =(doublecomplex *) sol;
    }

    // PrintMatrix(&A);
    CS           =atoi(argv[2]);
    LEVEL_OF_FILL=atoi(argv[3]);
    PA           =atoi(argv[4]);
    DROP_TOL     =atof(argv[5]);


    
    //-------------------------------------------------------------------------
    //------------------ set up user options for block ILU --------------------
    //-------------------------------------------------------------------------
    JanusDefaultOptions(&options);
    /* now the user may modify some of the options with respect
       to her (his) needs */
    /* turn matching on/off, default: on */
    options.matching=1;
    /* which reordering to take:
       PERM_NONE, PERM_MTMETIS, PERM_AMD, PERM_RCM
       default: PERM_MTMETIS */
    options.ordering=PERM_MTMETIS; // PERM_AMD; // PERM_NONE; // 
    /* drop tolerance, default: 1e-3 */
    options.droptol=DROP_TOL;
    /* cosine-based blocking turned on/off, default: off */
    options.cosine=CS;                 
    /* blocking turned on/off, default: BLOCK_ILU1T,
       alternatives: BLOCK_NONE|BLOCK_ILU1T|BLOCK_ILUPT|BLOCK_SUPERNODES */
    if (LEVEL_OF_FILL<0)
       options.blocking_strategy=BLOCK_SUPERNODES;
    else if (LEVEL_OF_FILL==0)
       options.blocking_strategy=BLOCK_NONE;
    else if (LEVEL_OF_FILL==1)
       options.blocking_strategy=BLOCK_ILU1T;
    else
       options.blocking_strategy=BLOCK_ILUPT;
    /* progressively aggreate blocks during the factorization
       on/off, default: on */
    options.progressive_aggregation=PA;
    /* allow diagonal block perturbations, default: on */
    options.perturbation=PERTURB; 
    /* work with symmetric permutations only and attempt to support 2x2 pivots
       default: off */
    options.symmetric_structure=0;
    // invert diagonal blocks (alternatively use factorization)
    options.invert_blocks=1;
    // level of fill, if ILUPT is selected
    options.level_of_fill=LEVEL_OF_FILL;
    //-------------------------------------------------------------------------
    //---------------- END set up user options for block ILU ------------------
    //-------------------------------------------------------------------------


        
    if (options.matching)
       printf("use maximum weight matching\n");
    else
       printf(" no maximum weight matching\n");
    if (options.cosine)
       printf("use cosine-based blocking\n");
    else
       printf(" no cosine-based blocking\n");
    if (options.ordering==PERM_MTMETIS)
       printf("use MT-METIS reordering\n");
    else if (options.ordering==PERM_METISN)
       printf("use METIS-NODE reordering\n");
    else if (options.ordering==PERM_METISE)
       printf("use METIS-EDGE reordering\n");
    else if (options.ordering==PERM_AMD)
       printf("use AMD reordering\n");
    else if (options.ordering==PERM_RCM)
       printf("use RCM reordering\n");
    else 
       printf("no reordering used\n");
    printf("use drop tolerance %8.1le\n",options.droptol);
    if (options.blocking_strategy==BLOCK_SUPERNODES) {
       printf("use supernodal blocking\n");
    }
    else if (options.blocking_strategy==BLOCK_ILU1T) {
       if (A.issymmetric||A.ishermitian)
	  if (A.isdefinite)
	     printf("use ICHOL(1,%8.1le) blocking\n",DROP_TOL);
	  else
	     printf("use BILDL(1,%8.1le) blocking\n",DROP_TOL);
       else
	  printf("use ILU(1,%8.1le)   blocking\n",DROP_TOL);
    }
    else if (options.blocking_strategy==BLOCK_ILUPT) {
       printf("use ILU(%2ld,%8.1le)  blocking\n",LEVEL_OF_FILL,DROP_TOL);
    }
    else // no blocking
       printf(" no blocking\n");
    if (options.progressive_aggregation)
       printf("use progressive_aggregation\n");
    else
       printf(" no progressive_aggregation\n");
    if (options.perturbation)
       printf("use perturbation\n");
    else
       printf(" no perturbation\n");
    if (!A.ishermitian && !A.issymmetric)
       if (options.symmetric_structure)
	  printf("use symmetric permutations and blocks\n");
       else
	  printf("non-symmetric permutations\n");
    if (options.invert_blocks)
       printf("use inverse diagonal blocks\n");
    else
       printf(" no inverse diagonal blocks\n");
    printf("\n");

    
#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
    //-------------------------------------------------------------------------
    //---------------------------- compute block ILU --------------------------
    //-------------------------------------------------------------------------
    JanusInit(&PREC);
    ierr=JanusFactor(&A, &PREC, options);
    //-------------------------------------------------------------------------
    //-------------------------- END compute block ILU ------------------------
    //-------------------------------------------------------------------------
#ifdef _PROFILING_
    time_bilu=omp_get_wtime()-timeBegin;
#endif


    
    if (ierr) {
       printf("computation of block incomplete factorization failed\n");

       JanusDelete(&PREC);
       SparseDelete(&A);
       free(sol);
       free(rhs);
       
       return (-1);
    }
    else {
       // printf("computation of block incomplete ILU successfully completed\n");
       //fflush(stdout);
    }

    if (PREC.isdefinite)
       printf("matrix is Hermitian positive definite\n");
    printf("log(det(A))=%12.4le+%12.4lei\n",PREC.logdet.r,PREC.logdet.i);
    
    // update number of blocks
    nblocks=PREC.BiD->nblocks;
    // some block statistics
    JanusNnz(&PREC, &j, &mx, &mu, &stddev);
    if (A.ishermitian || A.issymmetric) {
       // count true number of diagonal entries
       m=0;
       for (i=0; i<n; i++) 
	   for (k=0; k<A.ncol[i]; k++) {
	       l=A.rowind[i][k];
	       // diagonal entry exists
	       if (l==i) {
		  m++;
		  break;
	       } // end if
	   } // end for k
       // remember that for efficiency reasons only the lower triangular part of A is stored
       // JanusNnz acts in the symmetric case as if BUT=BL were present, though it is not stored
       fill_bilu=(double)j/(2*nz-m);
    }
    else
       fill_bilu=(double)j/nz;
    printf("block incomplete factorization, relative fill %8.1le\n",fill_bilu);
    printf("number of blocks: %ld, maximum block size: %ld, average %4.1lf(+-%8.1le)\n",nblocks,mx,mu,stddev);
    fflush(stdout);


    // artificial right hand side rhs=A*1 and initial guess 0
    if (A.isreal)
       for (i=0; i<n; i++) 
	   px[i]=1.0;	   
    else
       for (i=0; i<n; i++) {
	   pz[i].r=1.0;
	   pz[i].i=0.0;	   
       } // end for i
    MatVec(&A,sol,rhs);
    if (A.isreal)
       for (i=0; i<n; i++)
	   px[i]=0.0;	   
    else 
       for (i=0; i<n; i++) 
	   pz[i].r=pz[i].i=0.0;	   
   
    
#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
   
    //-------------------------------------------------------------------------
    // ---------- Krylov subspace solver preconditioned by BILU ---------------
    // ------------------------------------------------------------------------
    ierr=JanusSolver(&A,&PREC,rhs,sol,K_RESTART,RES_TOL,MAX_IT,&iter);
    // ------------------------------------------------------------------------
    // -------- END Krylov subspace solver preconditioned by BILU -------------
    // ------------------------------------------------------------------------
#ifdef _PROFILING_
    time_krylov=omp_get_wtime()-timeBegin;
#endif

    v=0.0;
    for (i=0; i<n; i++) {
        if (A.isreal)
	   v=MAX(fabs(px[i]-1.0),v);
	else
	   v=MAX(sqrt((pz[i].r-1.0)*(pz[i].r-1.0)+pz[i].i*pz[i].i),v);
    }
    printf("maximum error: %8.1le\n",v);fflush(stdout);


    
    // why did Krylov subspace solver stop?
    // error?
    if (ierr) {
       if (ierr==-1)
	  printf("too many iteration steps\n");
       else if (i==-2)
	  printf("not enough work space\n");
       else if (i==-3)
	  printf("algorithm breaks down\n");
       else 
	  printf("algorithm stops with error code %ld\n",ierr);
    }
    else { // success
       if (iter>=MAX_IT)
	  printf("preconditioned Krylov subspace solver stopped after %ld iteration steps\n", iter);
       else
	  printf("preconditioned Krylov subspace solver successfully completed with %ld iteration steps\n", iter);
    }

    
    
    //-------------------------------------------------------------------------
    // ------------------- clean up block triangular factors ------------------
    // ------------------------------------------------------------------------
    JanusDelete(&PREC);
    //-------------------------------------------------------------------------
    // ----------------- END clean up block triangular factors ----------------
    // ------------------------------------------------------------------------


    
#ifdef _PROFILING_
    time_total=omp_get_wtime()-time_total;
    printf("preconditioned iterative solver profiling summary\n");
    if (A.issymmetric||A.ishermitian)
       if (A.isdefinite) {
 	  printf(" (a) BICHOL                   %8.1le\n",time_bilu);
	  printf(" (b) CG       subspace solver %8.1le\n",time_krylov);
       }
       else {
	  printf(" (a) BILDL                    %8.1le\n",time_bilu);
	  printf(" (b) SQMR     subspace solver %8.1le\n",time_krylov);
       }
    else {
          printf(" (a) BILU                      %8.1le\n",time_bilu);
	  printf(" (b) GMRES(%2ld) subspace solver %8.1le\n",K_RESTART, time_krylov);
    }
    printf("Total solution time %8.1le\n\n",time_total);

    fflush(stdout);
#endif



    // release solution, right hand side and initial matrix
    free(sol);
    free(rhs);
    SparseDelete(&A);
       
    

    printf("-------------- %s                              ILU[sec]  rel.fill KRYLOV[sec] iter.steps  mxblk   avblk   stdv\n",argv[1]);
    printf("summary matrix %s,BILU(%1ld%1ld%1ld),droptol=%8.1le:  %8.1le  %8.1le   %8.1le  %5ld     %5ld %8.1lf  %8.1le\n",
	   argv[1],CS,LEVEL_OF_FILL,PA,DROP_TOL,
	   time_bilu,fill_bilu,time_krylov,iter,mx,mu,stddev);


    
    return (0);
} // end main

/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	August 29, 2017. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/


   $Id: janusdriver.c 6263 2020-05-15 07:04:30Z bolle $
*/
