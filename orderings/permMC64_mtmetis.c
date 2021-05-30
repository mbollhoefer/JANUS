#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include <blas.h>
#include <hsl.h>
#include <mtmetis.h>
#include <janus.h>
#include <ilupackmacros.h>

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

//#define PRINT_MEM
//#define PRINT_INFO
#define MEGA      1048576.0
#define IMEGA      262144.0
#define DMEGA      131072.0

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))


/*
 scaling and permutation driver, followed by a symmetric reordering.

 here 1) maximum weight matching (MC64)
      2) MeTis multilevel nested dissection by nodes


 Given an n x n matrix A in compressed sparse row format this routine computes
 row and column scalings as well as row and column permutations such that

 A -> Dr * A * Dc,
 
 where Dr and Dc are diagonal matrices stored in the vectors proscal and 
 pcoscal. The scaling is explicitly applied to A.

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
integer PERMMC64MTMETIS(CSRMAT A, FLOAT *proscal, FLOAT *pcoscal, integer *p, integer *invq, integer *nB, 
			ILUPACKPARAM *param)
{
/*
    A          n x n matrix in compressed sparse row format
               A is altered by the routine since the scalings are directly
               applied to A. Hopefully the scaling are only powers of 2,
               which is the case for all scaling routines from ILUPACK.

    proscal,   vectors of length n that store the row and column scalings Dr
    pcoscal    and Dc
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

        July/November 2003, November 2004, July 2005, April 2009. ILUPACK V2.3

    Notice:

	Copyright (c) 2009 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
 */

    integer  i,ii,j,k,l,m,imx,n=A.nr,liw,*iw,ldw,info[10],icntl[10],*cperm,job,
         iflag, num,nz,ierr=0, nnz,iwlen, scale, pow,
         flags, *ptr, *vwgt;
    REALS mxr, lg2=1.0/log((double)2.0), fpow, *b, *dw, val;
    FLOAT  *w,mx;
    CSRMAT B, C;
    integer bnd,bndt;
    double *mtmetis_options=NULL;
    size_t mem;
    REALS sm,sm2, sigma, sigmat, x;
     

    ierr=0;
    /* --------------------   END set up parameters   -------------------- */
    flags=param->ipar[6];

    // initial preprocessing
    if ((param->ipar[7]&(512+1024))==512) {
       scale=param->ipar[7]&(1+2+4);
    }
    // regular reordering
    else if ((param->ipar[7]&(512+1024))==1024) {
       scale=(param->ipar[7]&(8+16+32))>>3;
    }
    // final pivoting
    else if ((param->ipar[7]&(512+1024))==512+1024) {
       scale=(param->ipar[7]&(64+128+256))>>6;
    }
    else {
       scale=param->ipar[7]&(1+2+4);
    }
    /* if ((scale&4)==0)   => row scaling first
       else                   column scaling first
       
       if (scale&1)        => apply row scaling     \   ignored for 
       if (scale&2)        => apply column scaling  /   matching
    */


    // rescale the matrix a priori
    CLEAR(A.nc, pcoscal,1);
    // if desired apply row scaling and compute parameters for column scaling
    m=1;
    for (i=0; i<A.nc; i++)
    {
        j=A.ia[i]-1;
        k=A.ia[i+1]-1;
	// if row scaling should be applied before column scaling
	if ((scale&4)==0) {
	   l=k-j;
	   imx=I_AMAX(&l,A.a+j,&m);
	   mxr=FABS(A.a[j+imx-1]);

	   // compute nearest power of 2
	   fpow=log(mxr)*lg2;
	   if (fpow<0.0) {
	      pow=fpow-0.5;
	      fpow=1;
	      for (l=1; l<=-pow; l++)
		  fpow*=2.0;
	      mxr=fpow;
	   }
	   else {
	      pow=fpow+0.5;
	      fpow=1;
	      for (l=1; l<=pow; l++)
		  fpow*=2;
	      mxr=1.0/fpow;
	   }

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   proscal[i]=mx=mxr;
#else
	   mx.r=mxr;
	   mx.i=0.0;
	   proscal[i]=mx;
#endif
	   l=k-j;
	   SCAL(&l,&mx,A.a+j,&m);
	} // end if ((scale&4)==0)
        for (; j<k; j++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    if (pcoscal[A.ja[j]-1]<FABS(A.a[j])) {
	        pcoscal[A.ja[j]-1]=FABS(A.a[j]);
	    } // end if
#else
	    if (pcoscal[A.ja[j]-1].r<FABS(A.a[j])) {
	        pcoscal[A.ja[j]-1].r=FABS(A.a[j]);
	        pcoscal[A.ja[j]-1].i=0.0;
	    } // end if
#endif
	} // end for j
    } // end for i

    // invert column scaling parameters
    for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        mxr=pcoscal[i];
#else
        mxr=pcoscal[i].r;
#endif
	// compute nearest power of 2
	fpow=log(mxr)*lg2;
	if (fpow<0.0) {
	   pow=fpow-0.5;
	   // pow=floor(fpow);
	   fpow=1;
	   for (l=1; l<=-pow; l++)
	       fpow*=2.0;
	   mxr=fpow;
	}
	else {
	   pow=fpow+0.5;
	   // pow=ceil(fpow);
	   fpow=1;
	   for (l=1; l<=pow; l++)
	       fpow*=2;
	   mxr=1.0/fpow;
	}

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        pcoscal[i]=mxr;
#else
        pcoscal[i].r=mxr;
        pcoscal[i].i=0.0;
#endif
    } // end for i

    // apply column scaling
    for (i=0; i<n; i++)
    {
        j=A.ia[i]-1;
        k=A.ia[i+1]-1;
        for (; j<k; j++){
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    A.a[j]*=pcoscal[A.ja[j]-1];
#else
	    A.a[j].r*=pcoscal[A.ja[j]-1].r;
	    A.a[j].i*=pcoscal[A.ja[j]-1].r;
#endif
	} // end for j

	// if row scaling should be applied after column scaling
	if ((scale&4)!=0) {
	   j=A.ia[i]-1;
	   l=k-j;
	   imx=I_AMAX(&l,A.a+j,&m);
	   mxr=FABS(A.a[j+imx-1]);

	   // compute nearest power of 2
	   fpow=log(mxr)*lg2;
	   if (fpow<0.0) {
	      pow=fpow-0.5;
	      fpow=1;
	      for (l=1; l<=-pow; l++)
		  fpow*=2.0;
	      mxr=fpow;
	   }
	   else {
	      pow=fpow+0.5;
	      fpow=1;
	      for (l=1; l<=pow; l++)
		  fpow*=2;
	      mxr=1.0/fpow;
	   }

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   proscal[i]=mx=mxr;
#else
	   mx.r=mxr;
	   mx.i=0.0;
	   proscal[i]=mx;
#endif
	   l=k-j;
	   SCAL(&l,&mx,A.a+j,&m);
	} // end if ((scale&4)!=0) {
    }




    /* scale matrix such that the maximum entry in each column/row is one in
       absolute value */
    nz=A.ia[A.nc]-1;
    // cperm=(integer *) MALLOC((size_t)n*sizeof(integer),"permMC64_mtmetis:cperm");
    liw=5*n;
    if (liw<=param->niaux)
       iw=param->iaux;
    else
       iw=(integer *) MALLOC((size_t)liw*sizeof(integer),"permMC64_mtmetis:iw");
    ldw=3*n+nz;
    if (ldw<=param->ndaux)
       dw=(REALS *)param->daux;
    else
       dw=(REALS *)MALLOC((size_t)ldw*sizeof(FLOAT),"permMC64_mtmetis:dw");
    b=dw+3*n;

    /* -------------   construct permutation  ------------ */
#ifdef _DOUBLE_REAL_
    b=A.a;
#else
    for (i=0; i<nz; i++)
        b[i]=FABS(A.a[i]);
#endif
    MC64IR(icntl);
    job=5;
    MC64AR(&job,&n,&nz,A.ia,A.ja,b,&num,p,
	   &liw,iw,&ldw,dw,icntl,info);
    if (liw>param->niaux)
       iw=FREE(iw);

    for (i=0; i<n; i++) {
        mxr=exp(dw[i]);
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

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	dw[i]=mxr;
	pcoscal[i]*=dw[i];
#else
	dw[i]=mxr;
	pcoscal[i].r*=dw[i];
	pcoscal[i].i=0;
#endif
    } // end for i
    for (i=0; i<n; i++) {
        // wrong!
        // mxr=exp(dw[n+i])/dw[2*n+i];
        // correct:
        mxr=exp(dw[n+i]);

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

	dw[n+i]=mxr;
    } // end for i

       
    for (i=0; i<n; i++)	{
        j=A.ia[i]-1;
	l=A.ia[i+1]-1 - j;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	proscal[i]*=dw[n+i];
	SCAL(&l,dw+n+i,A.a+j,&m);
	l=l+j;
	for (; j<l; j++) {
	    A.a[j]*=dw[A.ja[j]-1];
	} // end for j
#else
	proscal[i].r*=dw[n+i];
	// SCAL(&l,dw+n+i,A.a+j,&m);
	l=l+j;
	for (; j<l; j++) {
	    A.a[j].r*=dw[n+i];
	    A.a[j].i*=dw[n+i];
	    A.a[j].r*=dw[A.ja[j]-1];
	    A.a[j].i*=dw[A.ja[j]-1];
	} // end for j
#endif
    }
    if (ldw>param->ndaux)
       dw=FREE(dw);







    if (flags&SYMMETRIC_STRUCTURE) {

       // compute symmetric weighted matching, i.e. break cycles into
       // 1x1 and 2x2 cycles
       // SYMMETRIC maximum weight matching interface
       if (3*A.nc+2>param->nibuff) {
	  param->nibuff=3*A.nc+2;
	  param->ibuff=(integer *)REALLOC(param->ibuff,param->nibuff*sizeof(integer),
					  "permMC64_mtmetis:param->ibuff");
       }
       if (3*A.nc>param->ndbuff) {
	  param->ndbuff=3*A.nc;
	  param->dbuff=(FLOAT *)REALLOC(param->dbuff,param->ndbuff*sizeof(FLOAT),
					"permMC64_mtmetis:param->dbuff");
       }
       l=GNLSYMMWM(A, p, param->ibuff, param->dbuff);


#ifdef PRINT_INFO
       printf("1x1: 1,...,%d\n",l);
       fflush(stdout);
       for (i=0; i<A.nr; i++) {
	   printf("%8d",p[i]);
       }
       printf("\n");
       fflush(stdout);
#endif


       // reorder the pattern of A by rows in order to apply a symbolic reordering

       // inverse permutation
       nz=A.ia[A.nr]-1;
       cperm=param->ibuff+2*A.nc+2;
       for (i=0; i<A.nc; i++)
	   cperm[p[i]-1]=i+1;
       C.nc=C.nr=A.nr;
       C.a=NULL;
       C.ia=param->ibuff+A.nc+1;
       if (param->niaux>=nz)
	  C.ja=param->iaux;
       else
	  C.ja=(integer *)MALLOC((size_t)nz*sizeof(integer), "permMC64_mtmetis:C.ja");
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
	       printf("%8d %8d %16.8le;...\n",i+1, C.ja[j-1],-1.0);
       printf("\n");
       fflush(stdout);
#endif


       // build symmetric undirected graph
       mem=(param->niaux>=k) ? param->niaux-k : 0;
       SETUPGRAPH(C,&B, param->ibuff, 
		  param->iaux+k,mem);
       // release memory for reordered pattern
       if (nz>param->niaux)
	  C.ja=FREE(C.ja);

       // check if MALLOC was needed to set up B.ja
       if (B.ia[B.nc]<0) {
	  iflag=0;
	  B.ia[B.nc]=-(B.ia[B.nc]);
       }
       else  {// nz+B.ia[B.nc]-1 entries from iaux were already in use 
	  iflag=B.ia[B.nc]-1; 
	  B.ja-=k;
	  // shift entries back by nz
	  for (i=0; i<iflag; i++) {
	      B.ja[i]=B.ja[i+k];
	  } // end for i
       }

       // get rid of the diagonal entries
       k=0;
       for (i=0; i<B.nc; i++) {
           m=k;
	   for (j=B.ia[i]-1; j<B.ia[i+1]-1;j++) {
	       // diagonal entry
	       if (B.ja[j]==i+1) {
		  // skip diagonal entry
		  k++;
	       }
	       else {
	          // shift diagonal entry
		  B.ja[j-k]=B.ja[j];
	       } // end if
	   } // end for j
	   B.ia[i]-=m;
       } // end for i
       B.ia[i]-=k;
	
       // reorder the system using multilevel nested dissection & minimum degree
       // switch from FORTRAN style to C-style
       for (i=0; i<=B.nc; i++)
	   B.ia[i]--;
       for (i=0; i<B.ia[B.nc]; i++)
	   B.ja[i]--;
       // default options
#ifndef _FULL_MTMETIS_
       int mt=omp_get_max_threads();
       //mexPrintf("max number of threads=%d\n",mt);fflush(stdout);
#ifdef _USE_MKL_
       //mexPrintf("MKL max number of threads=%d, dynamic=%d\n",mkl_get_max_threads(),mkl_get_dynamic());fflush(stdout);
       // mexPrintf("MKL domain max number of threads all =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_ALL));fflush(stdout);
       // mexPrintf("MKL domain max number of threads blas=%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_BLAS));fflush(stdout);
       // mexPrintf("MKL domain max number of threads vml =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_VML));fflush(stdout);
#endif
       omp_set_num_threads(MY_MTMETIS_THREADS);
#endif
       mtmetis_options=mtmetis_init_options();

       cperm=param->ibuff+A.nc+1;
       if (B.ia[B.nc]>1 && B.nc>16) {
	  vwgt=(integer *)malloc((size_t)B.nc*sizeof(integer));
	  for (i=0; i<B.nc; i++)
	      vwgt[i]=1;
	  // printf("Call MT-Metis after MC64\n");fflush(stdout);
	  MTMETIS_NodeND(&B.nc,B.ia,B.ja, vwgt,mtmetis_options,invq,cperm);
	  for (i=0; i<B.nc; i++) {
	      cperm[i]++;
	      invq[i]++;
	  } // end for i
	  free(vwgt);
       }
       else
	  for (i=0; i<B.nc; i++) 
	      cperm[i]=invq[i]=i+1;
       // switch from C-style to FORTRAN style 
       for (i=0; i<B.ia[B.nc]; i++)
	   B.ja[i]++;
       for (i=0; i<=B.nc; i++)
	   B.ia[i]++;
       if (!iflag)
	  B.ja=FREE(B.ja);
       free(mtmetis_options);
#ifndef _FULL_MTMETIS_
       omp_set_num_threads(mt);
       //mexPrintf("number of threads=%d\n",omp_get_num_threads());fflush(stdout);
       //mexPrintf("max number of threads=%d\n",omp_get_max_threads());fflush(stdout);
#ifdef _USE_MKL_
       // mexPrintf("MKL max number of threads=%d, dynamic=%d\n",mkl_get_max_threads(),mkl_get_dynamic());fflush(stdout);
       // mexPrintf("MKL domain max number of threads all =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_ALL));fflush(stdout);
       // mexPrintf("MKL domain max number of threads blas=%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_BLAS));fflush(stdout);
       // mexPrintf("MKL domain max number of threads vml =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_VML));fflush(stdout);
#endif
#endif
       /*
	 int MTMETIS_NodeND(
	 mtmetis_vtx_type const * nvtxs,
	 mtmetis_adj_type const * xadj,
	 mtmetis_vtx_type const * adjncy,
	 mtmetis_wgt_type const * vwgt,
	 double const * options,
	 mtmetis_pid_type * perm,
	 mtmetis_pid_type * iperm);
       */


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


    
       // build product permutation for the rows (maximum weighted matching 
       //                                         followed by
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


    }
    else { // generally structured matrix
       
       // reorder the pattern of A by rows in order to apply a symbolic reordering
       // allocate memory with respect to the specific permutation routine
       param->nibuff=MAX(param->nibuff,3*(size_t)A.nc+1);
       param->ibuff=(integer *) REALLOC(param->ibuff,param->nibuff*sizeof(integer),
					"permMC64_mtmetis:ibuff");

       // reorder the rows of A, at least the pattern of A
       C.nc=C.nr=A.nr;
       C.a=NULL;
       C.ia=param->ibuff;
       if (nz<=param->niaux)
	  C.ja=param->iaux;
       else
	  C.ja=(integer *)MALLOC((size_t)nz*sizeof(integer), "permMC64_mtmetis:C.ja");
       C.ia[0]=1;
       for (i=0; i<A.nc; i++) {
           ii=p[i]-1;
	   k=C.ia[i]-1;
	   for (j=A.ia[ii]; j<A.ia[ii+1]; j++) {
	       C.ja[k++]=A.ja[j-1];
	   }
	   C.ia[i+1]=C.ia[i]+ (A.ia[ii+1]-A.ia[ii]);
       }
       // build symmetric undirected graph
       // build B=|A|+|A^T|
       // param->ibuff+2n will also be used for B.ia!
       // remember that the leading 2n are already in use for symrcm
       // param->iaux will be used for B.ja if possible
       mem=(param->niaux>=C.ia[C.nc]-1) ? param->niaux-(C.ia[C.nc]-1) : 0;
       SETUPGRAPH(C,&B, param->ibuff+2*n, 
		  param->iaux+(C.ia[C.nc]-1),mem);

       // release memory for reordered pattern
       if (nz>param->niaux)
	  C.ja=FREE(C.ja);

       // check if MALLOC was needed to set up B.ja
       if (B.ia[B.nc]<0) {
          iflag=0;
	  B.ia[B.nc]=-(B.ia[B.nc]);
       }
       else  {// nz+B.ia[B.nc]-1 entries from iaux were already in use 
	  iflag=B.ia[B.nc]-1; 
       }

       // get rid of the diagonal entries
       k=0;
       for (i=0; i<A.nc; i++) {
           l=k;
	   for (j=B.ia[i]-1; j<B.ia[i+1]-1;j++) {
	       // diagonal entry
	       if (B.ja[j]==i+1) {
		  // skip diagonal entry
		 k++;
	       }
	       else {
		  // shift diagonal entry
		  B.ja[j-k]=B.ja[j];
	       } // end if
	   } // end for j
	   B.ia[i]-=l;
       } // end for i
       B.ia[i]-=k;
       
       // reorder the system using multilevel nested dissection & minimum degree
       // switch from FORTRAN to C-style
       for (i=0; i<=A.nc; i++)
	   B.ia[i]--;
       for (i=0; i<B.ia[A.nc]; i++)
	   B.ja[i]--;
       // default options
#ifndef _FULL_MTMETIS_
       int mt=omp_get_max_threads();
       // mexPrintf("max number of threads=%d\n",mt);fflush(stdout);
#ifdef _USE_MKL_
       // mexPrintf("MKL max number of threads=%d, dynamic=%d\n",mkl_get_max_threads(),mkl_get_dynamic());fflush(stdout);
       // mexPrintf("MKL domain max number of threads all =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_ALL));fflush(stdout);
       // mexPrintf("MKL domain max number of threads blas=%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_BLAS));fflush(stdout);
       // mexPrintf("MKL domain max number of threads vml =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_VML));fflush(stdout);
#endif
       omp_set_num_threads(MY_MTMETIS_THREADS);
#endif
       mtmetis_options=mtmetis_init_options();

       if (B.ia[A.nc]>1 && A.nc>16) {
	  vwgt=(integer *)malloc((size_t)A.nc*sizeof(integer));
	  for (i=0; i<A.nc; i++)
	      vwgt[i]=1;
	  // printf("Call MT-Metis after MC64\n");fflush(stdout);
	  MTMETIS_NodeND(&A.nc,B.ia,B.ja, vwgt,mtmetis_options,invq,param->ibuff+A.nc);
	  for (i=0; i<A.nc; i++) {
	      param->ibuff[A.nc+i]++;
	      invq[i]++;
	  } // end for i
	  free(vwgt);
       }
       else
	  for (i=0; i<A.nc; i++) 
	      param->ibuff[A.nc+i]=invq[i]=i+1;
       // switch from C-style to FORTRAN style 
       for (i=0; i<B.ia[A.nc]; i++)
	   B.ja[i]++;
       for (i=0; i<=A.nc; i++)
	   B.ia[i]++;
       free(mtmetis_options);
#ifndef _FULL_MTMETIS_
       omp_set_num_threads(mt);
       // mexPrintf("number of threads=%d\n",omp_get_num_threads());fflush(stdout);
       // mexPrintf("max number of threads=%d\n",omp_get_max_threads());fflush(stdout);
#ifdef _USE_MKL_
       // mexPrintf("MKL max number of threads=%d, dynamic=%d\n",mkl_get_max_threads(),mkl_get_dynamic());fflush(stdout);
       // mexPrintf("MKL domain max number of threads all =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_ALL));fflush(stdout);
       // mexPrintf("MKL domain max number of threads blas=%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_BLAS));fflush(stdout);
       // mexPrintf("MKL domain max number of threads vml =%d\n",mkl_domain_get_max_threads(MKL_DOMAIN_VML));fflush(stdout);
#endif
#endif
       /*
	 int MTMETIS_NodeND(
	 mtmetis_vtx_type const * nvtxs,
	 mtmetis_adj_type const * xadj,
	 mtmetis_vtx_type const * adjncy,
	 mtmetis_wgt_type const * vwgt,
	 double const * options,
	 mtmetis_pid_type * perm,
	 mtmetis_pid_type * iperm);
       */
       
       if (!iflag)
	  B.ja=FREE(B.ja);

       // build product permutation for the rows (minimum weighted matching followed 
       //                                         symmetric reordering)
       for (i=0; i<n; i++) {
	   param->ibuff[i]=p[param->ibuff[A.nc+i]-1];
       }
       // rewrite permutation
       for (i=0; i<n; i++) {
	   p[i]=param->ibuff[i];
       }
       
    } // end if-else (flags&SYMMETRIC_STRUCTURE)
    
    *nB=A.nc; 

    return (ierr);
} /* end permMC64_mtmetis */
