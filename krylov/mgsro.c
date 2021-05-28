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

#define IABS(A)         (((A)>=(0))?(A):(-(A)))
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout





void MGSRO(integer *full, integer *lda, integer *n, integer *m, integer *ind,
	   REALS *ops, FLOAT *vec, FLOAT *hh, integer *ierr)
{

  /*-----------------------------------------------------------------------
     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
               Orthogonalization
     The ind'th vector of VEC is orthogonalized against the rest of
     the vectors.

     The test for performing re-orthogonalization is performed for
     each indivadual vectors. If the cosine between the two vectors
     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
     performed. The norm of the 'new' vector is kept in variable NRM0,
     and updated after operating with each vector.

     full   -- .true. if it is necessary to orthogonalize the ind'th
               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
               .false. only orthogonalize againt vec(:,1:ind-1)
     lda    -- the leading dimension of VEC
     n      -- length of the vector in VEC
     m      -- number of vectors can be stored in VEC
     ind    -- index to the vector to be changed
     ops    -- operation counts
     vec    -- vector of LDA X M storing the vectors
     hh     -- coefficient of the orthogonalization
     ierr   -- error code
               0 : successful return
               -1: zero input vector
               -2: input vector contains abnormal numbers
               -3: input vector is a linear combination of others

     External routines used: FLOAT DISTDOT

     code taken from SPARSKIT of Yousef Saad.
     adapted by Matthias Bollhoefer for the complex case
  -----------------------------------------------------------------------
  */
  integer i,k,j,jj;
  REALS   nrm0, nrm1, thr, reorth;
  FLOAT   fct;

  // compute the norm of the input vector

  j=1;
  jj=1;
  nrm0=DISTDOT2(n,&vec[0+(*ind-1)**lda],&j);
  // write (6,'(A,1P,E12.4)') 'start Gram-Schmidt nrm0=',nrm0
  ops=ops+*n+*n;
  thr=nrm0*reorth;
  if (nrm0<=0.0) {
     *ierr=-1;
     return;
  }
  else if (nrm0>0.0 && 1.0/nrm0>0.0)
     *ierr=0;
  else {
     *ierr=-2;
     return;
  } // end if-else if-else

  // Modified Gram-Schmidt loop
  if (*full) {
     // write (6,'(A,I4,A,I4)')
     //     +        '1st Gram-Schmidt loop',ind+1,',...,',m
     for (i=*ind; i<*m; i++) {
         // write (6,'(A,I4)') ' Gram-Schmidt step',i+1
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
         fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
#ifdef _USE_MKL_
         DISTDOT(&fct,n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
         fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#endif
#endif
	 hh[i]=fct;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 fct=-fct;
#else
	 fct.r=-fct.r; fct.i=-fct.i;
#endif
	 AXPY(n,&fct,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 fct=-fct;
#else
	 fct.r=-fct.r; fct.i=-fct.i;
#endif
	 ops=ops+4**n+2;
	 if (FABS2(fct)>thr) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
#ifdef _USE_MKL_
	    DISTDOT(&fct,n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
	    fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#endif
#endif
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    hh[i]=hh[i]+fct;
	    fct=-fct;
#else
	    hh[i].r=hh[i].r+fct.r;
	    hh[i].i=hh[i].i+fct.i;
	    fct.r=-fct.r; fct.i=-fct.i;
#endif
	    AXPY(n,&fct,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    fct=-fct;
#else
	    fct.r=-fct.r; fct.i=-fct.i;
#endif
	    ops=ops+4**n+1;
	 } // end if
	 nrm0=nrm0-FABS2(hh[i]);
	 if (nrm0<0.0) nrm0=0.0;
	 thr=nrm0*reorth;
     } // end for i
  } // end if

  // write (6,'(A,I4,A,I4)') 
  //   +     '2nd Gram-Schmidt loop',1,',...,',ind-1
  for (i=0; i<*ind-1; i++) {
      // write (6,'(A,I4)') '2nd Gram-Schmidt step',i+1
      // write (6,'(1P,E12.4)') dznrm2(n,vec(1,ind),j)
      // write (6,'(1P,E12.4)') dznrm2(n,vec(1,i),j)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
#ifdef _USE_MKL_
      DISTDOT(&fct,n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
      fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#endif
#endif
      // write (6,'(A,1P,E12.4)') 'fct=',fct
      hh[i]=fct;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      fct=-fct;
#else
      fct.r=-fct.r; fct.i=-fct.i;
#endif
      AXPY(n,&fct,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
      // write (6,'(A)') 'axpy completed'
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      fct=-fct;
#else
      fct.r=-fct.r; fct.i=-fct.i;
#endif
      ops=ops+4**n+2;
      if (FABS2(fct)>thr) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
#ifdef _USE_MKL_
	 DISTDOT(&fct,n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#else
	 fct=DISTDOT( n,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
#endif
#endif
	 // write (6,'(A,1P,E12.4)') 'fct=',fct
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 hh[i]=hh[i]+fct;
	 fct=-fct;
#else
	 hh[i].r=hh[i].r+fct.r;
	 hh[i].i=hh[i].i+fct.i;
	 fct.r=-fct.r; fct.i=-fct.i;
#endif
	 AXPY(n,&fct,&vec[0+i**lda],&j,&vec[0+(*ind-1)**lda],&jj);
	 // write (6,'(A)') 'axpy completed'
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 fct=-fct;
#else
	 fct.r=-fct.r; fct.i=-fct.i;
#endif
	 ops=ops+4**n+1;
      } // end if
      nrm0=nrm0-FABS2(hh[i]);
      if (nrm0<0.0) nrm0=0.0;
      thr=nrm0*reorth;
  } // end for i

  // test the resulting vector

  nrm1=sqrt(DISTDOT2(n,&vec[0+(*ind-1)**lda],&j));
  // write (6,'(A,1P,E12.4)') ' norm resulting vector',nrm1
  ops=ops+*n+*n;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
  hh[*ind-1]=nrm1;
#else
  hh[*ind-1].r=nrm1; hh[*ind-1].i=0.0;
#endif
  if (nrm1<=0.0) {
     *ierr=-3;
     return;
  } // endif

  // scale the resulting vector
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
  fct=1.0/nrm1;
#else
  fct.r=1.0/nrm1; fct.i=0.0;
#endif
  SCAL(n,&fct,&vec[0+(*ind-1)**lda],&j);
  ops=ops+*n+1;

  // normal return

  *ierr=0;
  return;
} // end mgsro
