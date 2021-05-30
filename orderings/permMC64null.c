#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <hsl.h>
#include <janus.h>

#include <ilupackmacros.h>

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
integer PERMMC64NULL(CSRMAT A, FLOAT *proscal, FLOAT *pcoscal, integer *p, integer *invq, integer *nB, 
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
         iflag, num,nz,ierr=0, nnz,iwlen,ncmpa,options[10],numflag, scale, pow,
         flags, *ptr;
    REALS mxr, lg2=1.0/log((double)2.0), fpow, *b, *dw, val;
    FLOAT  *w,mx;
    CSRMAT B, C;
    integer bnd,bndt;
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
    // cperm=(integer *) MALLOC((size_t)n*sizeof(integer),"permMC64_null:cperm");
    liw=5*n;
    if (liw<=param->niaux)
       iw=param->iaux;
    else
       iw=(integer *) MALLOC((size_t)liw*sizeof(integer),"permMC64_null:iw");
    ldw=3*n+nz;
    if (ldw<=param->ndaux)
       dw=(REALS *)param->daux;
    else
       dw=(REALS *)MALLOC((size_t)ldw*sizeof(FLOAT),"permMC64_null:dw");
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
        // wrong! mxr=exp(dw[n+i])/dw[2*n+i];
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
       // compute symmetric weighted matching, i.e. break cycles into
       // 1x1 and 2x2 cycles
       // SYMMETRIC maximum weight matching interface
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



    
       // build product permutation for the rows (maximum weighted matching 
       //                                         followed by
       //                                         symmetric reordering)
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
       
       for (i=0; i<A.nc; i++) {
	   invq[i]=i+1;
       } // end for

       /*
       for (i=0; i<A.nc; i++) {
	   printf("%8d",invq[i]);
       } // end for
       printf("\n");
       fflush(stdout);
       */

    } // end if-else (flags&SYMMETRIC_STRUCTURE)
    
    *nB=A.nc; 

       /*
       for (i=0; i<A.nc; i++) {
	   printf("%8.1le",proscal[i]);
       } // end for
       printf("\n");
       for (i=0; i<A.nc; i++) {
	   printf("%8.1le",pcoscal[i]);
       } // end for
       printf("\n");
       fflush(stdout);
       */

    return (ierr);
} /* end permMC64_null */
