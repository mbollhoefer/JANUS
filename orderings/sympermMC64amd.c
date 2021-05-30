#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <hsl.h>
#include <amd.h>
#include <janus.h>

#include <ilupackmacros.h>

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout


// #define PRINT_INFO


#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)

#ifdef _SKEW_MATRIX_

#ifdef _DOUBLE_REAL_
#define MYSYMMWM             DSSMsmwm
#define MYSYMPERMMC64AMD     DSSMperm_mc64_amd

#else
#define MYSYMMWM             SSSMsmwm
#define MYSYMPERMMC64AMD     SSSMperm_mc64_amd
#endif

#define SKEW(A)      (-(A))

#else

#ifdef _DOUBLE_REAL_
#define MYSYMMWM             DSYMsmwm
#define MYSYMPERMMC64AMD     DSYMperm_mc64_amd
#else
#define MYSYMMWM             SSYMsmwm
#define MYSYMPERMMC64AMD     SSYMperm_mc64_amd
#endif

#define SKEW(A)      (A)
#endif
// end _SKEW_MATRIX_


#else

#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)

#ifdef _SKEW_MATRIX_
#define SKEW(A)      (-(A))

#ifdef _SINGLE_COMPLEX_
#define MYSYMMWM             CSSMsmwm
#define MYSYMPERMMC64AMD     CSSMperm_mc64_amd
#else
#define MYSYMMWM             ZSSMsmwm
#define MYSYMPERMMC64AMD     ZSSMperm_mc64_amd
#endif

#else
#define SKEW(A)      (A)

#ifdef _SINGLE_COMPLEX_
#define MYSYMMWM            CSYMsmwm
#define MYSYMPERMMC64AMD    CSYMperm_mc64_amd
#else
#define MYSYMMWM            ZSYMsmwm
#define MYSYMPERMMC64AMD    ZSYMperm_mc64_amd
#endif

#endif
// end _SKEW_MATRIX_


#else
#define CONJG(A)     (-(A))

#ifdef _SKEW_MATRIX_
#define SKEW(A)      (-(A))

#ifdef _SINGLE_COMPLEX_
#define MYSYMMWM            CSHRsmwm
#define MYSYMPERMMC64AMD    CSHRperm_mc64_amd
#else
#define MYSYMMWM            ZSHRsmwm
#define MYSYMPERMMC64AMD    ZSHRperm_mc64_amd
#endif

#else
#define SKEW(A)      (A)

#ifdef _SINGLE_COMPLEX_
#define MYSYMMWM            CHERsmwm
#define MYSYMPERMMC64AMD    CHERperm_mc64_amd
#else
#define MYSYMMWM            ZHERsmwm
#define MYSYMPERMMC64AMD    ZHERperm_mc64_amd
#endif

#endif
// end _SKEW_MATRIX_

#endif
// end _COMPLEX_SYMMETRIC_

#endif
// end _DOUBLE_REAL_


/*
 scaling and permutation driver, followed by a symmetric reordering.

 here 1) maximum weight matching (MC64)
      3) associated symmetric block partitioning and block reordering by
         I.S. Duff and S. Pralet
      2) AMD (approximate minimum degree) from UMFPACK V4.3


 Given an n x n matrix A in compressed sparse row format this routine computes
 row and column scalings as well as row and column permutations such that

 A -> Dr * A * Dc,
 
 where Dr and Dc are diagonal matrices stored in the vectors proscal and 
 pcolscale. The scaling is explicitly applied to A.

 The permutation p,invq refer to permutation matrices P and Q^{-1} such  that

  P^T*(Dr * A * Dc)*Q 

 is hopefully better suited for the ILU.

 The routine returns a leading block size nB<=n
                        /B F\ 
  P^T*(Dr * A * Dc)*Q = |   | 
                        \E C/
 In this case the permutation recommends that the ILU is only applied to the
 leading block B of size nB x nB
 */
integer MYSYMPERMMC64AMD(CSRMAT A, FLOAT *prowscale, FLOAT *pcolscale, 
		     integer *p, integer *invq, integer *nB, ILUPACKPARAM *param)
{
/*
    A          n x n matrix in compressed sparse row format
               A is altered by the routine since the scalings are directly
               applied to A. Hopefully the scaling are only powers of 2,
               which is the case for all scaling routines from ILUPACK.

    prowscale,   vectors of length n that store the row and column scalings Dr
    pcolscale    and Dc
               A -> Dr*A*Dc

    p,invq     permutation vectors p and q^{-1} of length n which refer to
               row / column permutation matrices P and Q^{-1}.
               P^T*(Dr * A * Dc)*Q hopefully better suited for ILU

    nB         leading blocksize nB
                                     /B F\ 
               P^T*(Dr * A * Dc)*Q = |   |,  B matrix of size nB x nB
                                     \E C/
               nB is the recommended block size for the application of an ILU

    param      ILUPACK parameters


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        November 2004, May 2005, April 2009. ILUPACK V2.3

    Notice:

	Copyright (c) 2009 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
 */ 

    integer  i,ii,j,k,l,m,imx,n=A.nr,liw,*iw,ldw,info[10],icntl[10],*cperm,job,
         iflag,
         num,nz,ierr=0,options[10], numflag, scale, pow,nrm, *ia, *ja,bnd,bndt;
    REALS mxr, lg2=1.0/log((double)2.0), fpow, *dw, a, *b,
          sm,sm2, sigma, sigmat, x;
    FLOAT  *w,mx;
    double Control[AMD_CONTROL], Info[AMD_INFO];
    size_t mem;
    CSRMAT B, C;

    // mexPrintf("entering SYMPERMMC64AMD\n");fflush(stdout);
    ierr=0;
    /* --------------------   END set up parameters   -------------------- */


    param->ndbuff=MAX(param->ndbuff,2*(size_t)A.nc);
    param->dbuff=(FLOAT *) REALLOC(param->dbuff,param->ndbuff*sizeof(FLOAT),
				   "sympermMC64_amd:dbuff");

#include "scaleprefix.c"

    // use unsymmetric version of A for maximum weight matching
    // build B=|A|+|A^T|
    param->nibuff=MAX(param->nibuff,3*(size_t)A.nc+2);
    param->ibuff=(integer *) REALLOC(param->ibuff,param->nibuff*sizeof(integer),
				 "sympermMC64_amd:ibuff");
    // param->ibuff (1,...,2n+1) will be used indicate the pruned parts of A, but 
    // also be used for B.ia!
    // the leading 5n entries of param->iaux will be used for MC64, if possible
    // dbuff (n) will be used to compute the maximum absolute entry per row.
    mem=(param->niaux>=5*n) ? param->niaux-5*n : 0;
    SETUPGRAPH_EPSILON(A,&B,(REALS)MC64THRESHOLD,
		       (REALS *)param->dbuff, param->ibuff, 
		       param->iaux+5*n, mem);

    // check if MALLOC was needed to set up B.ja
    if (B.ia[B.nc]<0) {
       iflag=0;
       B.ia[B.nc]=-(B.ia[B.nc]);
       // printf("%ld,%u, memory for B.ja allocated\n",(B.ia[B.nc]-1),mem); fflush(stdout);
    }
    else  {// nz entries from iaux+5n are already in use 
       iflag=B.ia[B.nc]-1; 
       // printf("%ld,%u, B.ja remapped\n",(B.ia[B.nc]-1),mem); fflush(stdout);
    }

    // number of nonzeros
    nz=B.ia[B.nc]-1;

    // printf("nnz(A)=%10ld, nnz(B)=%10ld\n",A.ia[A.nc]-1,B.ia[B.nc]-1);fflush(stdout);

    // memory for the unsymmetric version of A
    ldw=3*B.nc+nz;
    if (param->ndaux>=ldw) {
       dw=(REALS *)param->daux;
       // printf("%ld,%ld,remap dw\n",ldw,param->ndaux); fflush(stdout);
    }
    else {
       dw=(REALS *) MALLOC((size_t)ldw*sizeof(REALS),"sympermMC64amd:dw");
       // printf("%ld,%ld, allocate memory for dw\n",ldw,param->ndaux); fflush(stdout);
    }
    b=dw+3*B.nc;

    /*--------------------  actual copying  */
    iw=B.ia;
    for (i=0; i<A.nc; i++) {
        // copy indices
        for (j=A.ia[i]-1; j<param->ibuff[i]-1; j++) {
	    // current column index of row i in FORTRAN style
	    k=A.ja[j]-1;
	    a=FABS(A.a[j]);
	    b[iw[i]-1]=a;
	    // advance pointer
	    iw[i]++;
	    // avoid duplicate entries
	    if (k!=i) {
	       b[iw[k]-1]=a;
	       // advance pointer
	       iw[k]++;
	    } // end if
	} // end for j
    } // end for i
    // shift iw back
    for (i=A.nc; i>0; i--) 
        iw[i]=iw[i-1];
    iw[0]=1;


    // MC64 maximum weight matching interface
    liw=5*B.nc;
    if (param->niaux>=liw) {
       iw=param->iaux;
       // printf("%ld,%ld,remap iw\n",liw,param->niaux); fflush(stdout);
    }
    else {
       iw=(integer *) MALLOC((size_t)liw*sizeof(integer),"sympermMC64amd:iw");
       // printf("%ld,%ld,allocate memory for iw\n",liw,param->niaux); fflush(stdout);
    }

    /* -------------   construct permutation  ------------ */
    MC64IR(icntl);
    job=5;
    MC64AR(&job,&(B.nc),&nz,B.ia,B.ja,b,&num,p,
	   &liw,iw,&ldw,dw,icntl,info);

    if (param->niaux<liw) {
       iw=FREE(iw);
       // printf("release iw\n"); fflush(stdout);
    }
    else {
       // printf("%4d,%4d,%4d\n",iw[0],iw[1],iw[2]);fflush(stdout);
       // printf("iw mapped to parts of iaux\n"); fflush(stdout);
    }
    nz=A.ia[A.nc]-1;
    // values of the unsymmetric version are not longer needed
    if (!iflag) {
       B.ja=FREE(B.ja);
       // printf("release B.ja\n"); fflush(stdout);
    }
   
    // build left and right diagonal scaling derived from the matching
    for (i=0; i<n; i++) {

	mxr=sqrt(exp(dw[i]+dw[n+i]));
	// compute nearest power of 2
	fpow=log(mxr)*lg2;
	if (fpow<0.0) {
	   pow=fpow-0.5;
	   fpow=1;
	   for (l=1; l<=-pow; l++)
	       fpow*=2.0;
	   mxr=1.0/fpow;
	}
	else {
	   pow=fpow+0.5;
	   fpow=1;
	   for (l=1; l<=pow; l++)
	       fpow*=2;
	   mxr=fpow;
	}
	dw[i]=mxr;
    } // end for i

    // transfer scaling to the official interface vectors
    for (i=0; i<A.nr; i++){
#if defined _DOUBLE_REAL_ || defined  _SINGLE_REAL_
	pcolscale[i]*=dw[i];
#else
	pcolscale[i].r*=dw[i];
	pcolscale[i].i=0;
#endif
    } // end for i


    // rescale matrix explicitly
    for (i=0; i<n; i++){
        for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++){
#if defined _DOUBLE_REAL_ || defined  _SINGLE_REAL_
	    A.a[j]=dw[i]*A.a[j]*dw[A.ja[j]-1];
#else
	    A.a[j].r=dw[i]*A.a[j].r*dw[A.ja[j]-1];
	    A.a[j].i=dw[i]*A.a[j].i*dw[A.ja[j]-1];
#endif
	} // end for j
    } // end for i
    if (param->ndaux<ldw) {
       dw=FREE(dw);
       //printf("release dw\n"); fflush(stdout);
    }

#ifdef PRINT_INFO
    printf("matrix after rescaling\n");
    for (i=0; i<A.nc; i++)
        for (j=A.ia[i];j<A.ia[i+1];j++)
	    printf("%8d%8d%16.8le\n",i+1, A.ja[j-1],A.a[j-1]);
    printf("\n");
    fflush(stdout);
    for (i=0; i<A.nc; i++) 
        printf("%8d",p[i]);
    printf("\n");
    fflush(stdout);
#endif



    // if the problem is already block structured then we do not
    // break up the block structure by matchings but simply 
    // reorder the block structured and scaled system.
    if (param->flags&BLOCK_STRUCTURE) {
       for (i=0; i<=A.nc; i++) 
	   invq[i]=0;
       // count number of blocks
       k=0;
       for (i=0; i<A.nc; i++) {
	   j=param->indicator[i];
	   // maximum block
	   if (j>k)
	      k=j;
	   // increase number of entries of block j
	   invq[j]++;
       }
#ifdef PRINT_INFO
       printf("number of blocks %ld, block size:\n",k);
       for (i=0; i<=k; i++)
           printf("%8d",invq[i]);
       printf("\n\n");
       fflush(stdout);
#endif

       // switch to pointer structure
       for (i=0; i<k; i++) {
	   if (invq[i+1]<=0) {
	      printf("sympermMC64_amd: block structure is incorrect for block %ld\n",i+1);
	      exit(0);
	   }
	   invq[i+1]+=invq[i];
       }
#ifdef PRINT_INFO
       printf("number of blocks %ld, pointer structure:\n",k);
       for (i=0; i<=k; i++)
           printf("%8d",invq[i]);
       printf("\n\n");
       fflush(stdout);
#endif

       // rearrange matrix such that blocks are grouped together
       for (i=0; i<A.nc; i++) {
	   j=param->indicator[i];
	   p[invq[j-1]++]=i+1;	     
       }
#ifdef PRINT_INFO
       printf("rearranged entries\n");
       for (i=0; i<A.nc; i++)
           printf("%8d",p[i]);
       printf("\n\n");
       fflush(stdout);
#endif

       /*
       // shift pointer array forward
       for (i=k-1; i>0; i--)
	   invq[i+1]=invq[i];
       invq[0]=0;
       */

       cperm=param->ibuff+A.nc+1;
       // compute inverse permutation
       for (i=0; i<A.nc; i++)
	   cperm[p[i]-1]=i+1;

       // compute reordered matrix
       C.nc=C.nr=A.nr;
       C.a=NULL;
       C.ia=param->ibuff;
       if (param->niaux>=nz) {
	  C.ja=param->iaux;
	  // printf("C.ja taken from iaux\n");
       }
       else {
	  C.ja=(integer *)MALLOC((size_t)nz*sizeof(integer), "sympermMC64_amd:C.ja");
	  // printf("C.ja allocated\n");
       }
       C.ia[0]=1;
       for (i=0; i<A.nc; i++) {
	   ii=p[i]-1;
	   k=C.ia[i]-1;
	   for (j=A.ia[ii]; j<A.ia[ii+1]; j++) {
	       C.ja[k++]=cperm[A.ja[j-1]-1];
	   }
	   C.ia[i+1]=C.ia[i]+ (A.ia[ii+1]-A.ia[ii]);
       }



       // transfer block structure to cperm
       for (i=0; i<A.nc; i++)
	   cperm[i]=param->indicator[p[i]-1];

#ifdef PRINT_INFO
       printf("block structure\n");
       for (i=0; i<A.nc; i++)
           printf("%8d",cperm[i]);
       printf("\n\n");
       fflush(stdout);
#endif

       // overwrite column indices
       for (i=0; i<nz; i++)
	   C.ja[i]=cperm[C.ja[i]-1];

       // set up buffer for flags
       for (i=0; i<A.nc; i++)
	   cperm[i]=0;
       
       // merge duplicate entries and rows 
       m=0;
       k=0;
       for (i=0; i<A.nc; ) {
	   // compute current block size
	   iflag=-1;
	   l=i+1;
	   while (iflag) {
	         if (l>=A.nc)
		    iflag=0;
		 else {
		    if (param->indicator[p[l]-1]!=param->indicator[p[i]-1])
		       iflag=0;
		    else
		       l++;
		 }
	   }
	   l=l-i;
	   // printf("i=%ld, l=%ld\n",i,l);

           // old start of row i
           j=C.ia[i]-1;
	   // new start of row i
	   C.ia[m]=k+1;
	   for (; j<C.ia[i+l]-1; j++) {
	       ii=C.ja[j]-1;
	       // entry is not duplicate yet
	       if (!cperm[ii]) {
		  // check mark
		  cperm[ii]=-1;
		  // shift entry
		  C.ja[k++]=ii+1;
	       }
	   }
	   // clear check mark array
	   for (j=C.ia[m]-1; j<k; j++) {
	       ii=C.ja[j]-1;
	       cperm[ii]=0;
	   }
	   i+=l;
	   m++;
       }
       C.ia[m]=k+1;
       C.nc=C.nr=m; 

#ifdef PRINT_INFO
       for (i=0; i<C.nc; i++)
	   for (j=C.ia[i];j<C.ia[i+1];j++)
	       printf("%8d%8d%16.8le\n",i+1, C.ja[j-1],-1.0);
       printf("\n");
       fflush(stdout);
#endif


       // setup matrix for AMD
       for (i=0; i<C.nc; i++) {
	   j=C.ia[i]-1;
	   k=C.ia[i+1]-1-j;
	   QQSORTI(C.ja+j,param->ibuff+A.nc+1,&k);
       }


       // convert from FORTRAN style to C-style
       for (i=0; i<=C.nc; i++)
	   C.ia[i]--;
       for (i=0; i<nz; i++)
	   C.ja[i]--;

       // reorder the system using approximate minimum degree
       cperm=invq;
#ifdef _LONG_INTEGER_
       amd_l_defaults(Control);
       ierr=amd_l_order(C.nc,C.ia,C.ja,cperm,Control,Info);
#else
       amd_defaults(Control);
       ierr=amd_order(C.nc,C.ia,C.ja,cperm,Control,Info);
#endif

       // printf("AMD symperm block ierr=%ld\n",ierr); fflush(stdout);
       if (param->niaux<nz) {
	  C.ja=FREE(C.ja);
	  //printf("release C.ja\n");fflush(stdout);
       }

       // copy cperm to a different buffer
       cperm=param->ibuff+A.nc+1;
       for (i=0; i<C.nc; i++) {
	   cperm[i]=++invq[i];
       } // end for



#ifdef PRINT_INFO
       printf("AMD block permutation\n");
       for (i=0; i<C.nc; i++)
	   printf("%8d",cperm[i]);
       printf("\n\n");
       fflush(stdout);
#endif



       // set up block pointer list
       for (i=0; i<=A.nc; i++) 
	   invq[i]=0;
       // count number of blocks
       k=0;
       for (i=0; i<A.nc; i++) {
	   j=param->indicator[i];
	   // maximum block
	   if (j>k)
	      k=j;
	   // increase number of entries of block j
	   invq[j]++;
       }
       // switch to pointer structure
       for (i=0; i<k; i++) {
	   invq[i+1]+=invq[i];
       }
       
       // prolongate permutation to its original size
       j=0;
       // l=cperm[0];
       // m=cperm[0];
       for (i=0;i<B.nc; i++)
	   for (k=invq[cperm[i]-1]; k<invq[cperm[i]]; k++) {
	       // if (cperm[i]-1<0) { printf("index %ld,%ld, index access below lower bound\n",i+1,cperm[i]); fflush(stdout);}
	       // if (cperm[i]>=A.nc) { printf("index %ld,%ld, index access above upper bound\n",i+1,cperm[i]+1); fflush(stdout); }
	       // l=MAX(l,cperm[i]);
	       // m=MIN(m,cperm[i]);
	       param->ibuff[j++]=k+1;
	   }

       // printf("range 0<= %ld,%ld<=%ld\n",l,m,B.nc); fflush(stdout);

#ifdef PRINT_INFO
       printf("prolongated block permutation\n");fflush(stdout);
       for (i=0; i<A.nc; i++)
	   printf("%8d",param->ibuff[i]);
       printf("\n\n");
       fflush(stdout);
#endif


    
       // build product permutation for the rows (maximum 
       // weighted matching followed symmetric reordering)
       for (i=0; i<n; i++) {
	   cperm[i]=p[param->ibuff[i]-1];
       }
       // rewrite permutation
       for (i=0; i<n; i++) {
	   p[i]=cperm[i];
       }
       // inverse permutation
       for (i=0; i<n; i++) {
	   invq[p[i]-1]=i+1;
       }



#ifdef PRINT_INFO
       printf("final scaling\n");fflush(stdout);
       for (i=0; i<A.nr; i++){
#if defined _DOUBLE_REAL_ || defined  _SINGLE_REAL_
	   printf("%16.8le", pcolscale[i]);
#else
	   printf("%16.8le", pcolscale[i].r);
#endif
       } // end for i
       printf("\n\n");
       fflush(stdout);
       printf("final permutation\n");fflush(stdout);
       for (i=0; i<A.nc; i++)
	   printf("%8d",p[i]);
       printf("\n\n");
       fflush(stdout);
#endif

       *nB=A.nc; 

       return (ierr);
    } // end BLOCK_STRUCTURE



    // compute symmetric weighted matching, i.e. break cycles into
    // 1x1 and 2x2 cycles
    l=MYSYMMWM(A, p, param->ibuff, param->dbuff);
    


#ifdef PRINT_INFO
    for (i=0; i<A.nr; i++)
      printf("%8d",p[i]);
    printf("\n");
    fflush(stdout);
#endif

 
    // reorder the rows and columns of A, at least the pattern of it
    // inverse permutation
    // ibuff (2n+1) entries are used for C.ia and cperm
    cperm=param->ibuff+A.nc+1;
    for (i=0; i<A.nc; i++)
        cperm[p[i]-1]=i+1;
    C.nc=C.nr=A.nr;
    C.a=NULL;
    C.ia=param->ibuff;
    if (param->niaux>=nz) {
       C.ja=param->iaux;
       // printf("%ld,%ld, remap C.ja\n",nz,param->niaux); fflush(stdout);
    }
    else {
       C.ja=(integer *)MALLOC((size_t)nz*sizeof(integer), "sympermMC64_amd:C.ja");
       // printf("%ld,%ld, allocate memory for C.ja\n",nz,param->niaux); fflush(stdout);
    }
    C.ia[0]=1;
    for (i=0; i<A.nc; i++) {
        ii=p[i]-1;
	k=C.ia[i]-1;
	for (j=A.ia[ii]; j<A.ia[ii+1]; j++) {
	    C.ja[k++]=cperm[A.ja[j-1]-1];
	}
	C.ia[i+1]=C.ia[i]+ (A.ia[ii+1]-A.ia[ii]);
    }


    // compress matrix
    // the leading l rows correspond to 1x1 pivots
    for (i=0; i<l; i++)
        cperm[i]=i+1;
    j=l;
    for (; i<A.nc; i+=2)
        cperm[i]=cperm[i+1]=++j;


    // shift column indices
    for (i=0; i<nz; i++)
        C.ja[i]=cperm[C.ja[i]-1];

    // buffer for flags
    for (i=0; i<A.nc; i++)
        cperm[i]=0;
    // skip duplicate entries in the leading l rows
    k=0;
    for (i=0; i<l; i++) {
        // old start of row i
        j=C.ia[i]-1;
	// new start of row i
	C.ia[i]=k+1;
        for (; j<C.ia[i+1]-1; j++) {
	    ii=C.ja[j]-1;
	    // entry is not duplicate yet
	    if (!cperm[ii]) {
	       // check mark
	       cperm[ii]=-1;
	       // shift entry
	       C.ja[k++]=ii+1;
	    }
	}
	// clear check mark array
        for (j=C.ia[i]-1; j<k; j++) {
	    ii=C.ja[j]-1;
	    cperm[ii]=0;
	}
    }
    // merge duplicate entries and rows in the remaining  n-l rows
    m=l;
    for (; i<A.nc; i+=2,m++) {
        // old start of row i
        j=C.ia[i]-1;
	// new start of row i
	C.ia[m]=k+1;
        for (; j<C.ia[i+2]-1; j++) {
	    ii=C.ja[j]-1;
	    // entry is not duplicate yet
	    if (!cperm[ii]) {
	       // check mark
	       cperm[ii]=-1;
	       // shift entry
	       C.ja[k++]=ii+1;
	    }
	}
	// clear check mark array
        for (j=C.ia[m]-1; j<k; j++) {
	    ii=C.ja[j]-1;
	    cperm[ii]=0;
	}
    }
    C.ia[m]=k+1;
    C.nc=C.nr=m;



#ifdef PRINT_INFO
    for (i=0; i<C.nc; i++)
        for (j=C.ia[i];j<C.ia[i+1];j++)
	    printf("%8d%8d%16.8le\n",i+1, C.ja[j-1],-1.0);
    printf("\n");
    fflush(stdout);
#endif



    // setup matrix for AMD
    for (i=0; i<C.nc; i++) {
        j=C.ia[i]-1;
	k=C.ia[i+1]-1-j;
	QQSORTI(C.ja+j,param->ibuff+A.nc+1,&k);
    }


    // convert from FORTRAN style to C-style
    for (i=0; i<=C.nc; i++)
        C.ia[i]--;
    for (i=0; i<nz; i++)
        C.ja[i]--;



    // reorder the system using approximate minimum degree
    cperm=invq;
#ifdef _LONG_INTEGER_
    amd_l_defaults(Control);
    ierr=amd_l_order(C.nc,C.ia,C.ja,cperm,Control,Info);
#else
    amd_defaults(Control);
    ierr=amd_order(C.nc,C.ia,C.ja,cperm,Control,Info);
#endif

    // printf("AMD symperm dynamic ierr=%ld\n",ierr); fflush(stdout);
    if (param->niaux<nz) {
       C.ja=FREE(C.ja);
       //printf("release C.ja\n");fflush(stdout);
    }

    for (i=0; i<C.nc; i++) {
        cperm[i]++;
    } // end for




#ifdef PRINT_INFO
    for (i=0; i<C.nc; i++)
        printf("%8d",cperm[i]);
    printf("\n\n");
    fflush(stdout);
#endif



    // prolongate permutation to its original size
    j=0;
    for (i=0;i<C.nc; i++)
        // 1x1 block
        if (cperm[i]<=l)
	   param->ibuff[j++]=cperm[i];
        else { // 2x2 block
	   param->ibuff[j]=2*cperm[i]-l-1;
	   param->ibuff[j+1]=param->ibuff[j]+1;
	   j+=2;
	}



#ifdef PRINT_INFO
    printf("block permutation\n");fflush(stdout);
    for (i=0; i<A.nc; i++)
        printf("%8d",param->ibuff[i]);
    printf("\n\n");
    fflush(stdout);
#endif


    
    // build product permutation for the rows (minimum weighted matching followed 
    //                                         symmetric reordering)
    for (i=0; i<n; i++) {
        cperm[i]=p[param->ibuff[i]-1];
    }
    // rewrite permutation
    for (i=0; i<n; i++) {
	p[i]=cperm[i];
    }
    // inverse permutation
    for (i=0; i<n; i++) {
	invq[p[i]-1]=i+1;
    }



#ifdef PRINT_INFO
    printf("final permutation\n");fflush(stdout);
    for (i=0; i<A.nc; i++)
        printf("%8d",p[i]);
    printf("\n\n");
    fflush(stdout);
#endif

    *nB=A.nc; 

    return (ierr);
} /* end sympermMC64_amd */
