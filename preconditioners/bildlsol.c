/* $Id: bildlsol.c 4169 2018-04-26 16:34:10Z bolle $ */
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

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)
#define TRANSA        "t"

#ifdef _DOUBLE_REAL_
#define MYBILDLSOL DSYMbilusol
#define MYSYTRS    dsytrs
#else
#define MYBILDLSOL SSYMbilusol
#define MYSYTRS    ssytrs
#endif
// endif-else _DOUBLE_REAL_

#else
// !defined _DOUBLE_REAL_ and !defined _SINGLE_REAL_

#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)
#define TRANSA        "t"

#ifdef _SINGLE_COMPLEX_
#define MYBILDLSOL CSYMbilusol
#define MYSYTRS    csytrs
#else
#define MYBILDLSOL ZSYMbilusol
#define MYSYTRS    zsytrs
#endif
// end if-else _SINGLE_COMPLEX_

#else
// !defined _COMPLEX_SYMMETRIC_

#define CONJG(A)     (-(A))
#define TRANSA        "c"

#ifdef _SINGLE_COMPLEX_
#define MYBILDLSOL CHERbilusol
#define MYSYTRS    chetrs
#else
#define MYBILDLSOL ZHERbilusol
#define MYSYTRS    zhetrs
#endif
// end if-else _SINGLE_COMPLEX_


#endif
// end if-else !defined _COMPLEX_SYMMETRIC_

#endif
// end if-else defined _DOUBLE_REAL_ or defined _SINGLE_REAL_

// #define PRINT_INFO
// #define printf mexPrintf

void MYBILDLSOL(SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD, SPARSEBLOCKMATRIX *BUT, integer *pivots,
		FLOAT *rhssol, FLOAT *dbuff, integer n, integer m)
{  
  /*
    solve linear system(s) with block incomplete LDL^T factorization

    Input
    -----
    BL          lower block triangular matrix L with unit diagonal
    BiD         block diagonal matrix D^{-1} (LAPACK-factorized and inverted)
    BUT         not referenced
    rhssol      right hand side(s) on input
    n           number of rows
    m           number of columns

    Output
    ------
    rhssol      computed solution on output
  */

  integer i,j,k,l,r,tt,i2, // counters
          *pi;          // index array pointer
  
  FLOAT   *pL,*pD,*pb, // pointers to lower/diagonal/upper/buff
          *buff,       // buffer for the current block column/row of BL/BUT
          *psol,       // pointer to solution
          val,         // auxiliary value
          alpha, beta; // parameters for level-2/3-BLAS

  char    *transa, *transb; // strings for level-2/3-BLAS
  REALS   rval;
  integer invert_blocks=1, spd=0,ierr;


  // if pivots=0, then diagonal blocks are inverse matrices with which we could
  // simply multiply. Otherwise pivoting is required. There is one exception,
  // which is the SPD case which never requires pivoting. In order to distinguish
  // the inversion case from the factorized situation, we pass a dummy pivot s.t.
  // its leading entries is -(n+1). This refers in the SPD case to the Cholesky case
  if (pivots!=NULL) {
     // a pivoting vector is passed -> LDL^T=P^TAP or LL^T=A
     invert_blocks=0;
     // LL^T=A, Cholesky case
     if (pivots[0]==-(n+1))
        spd=1;
  } // end if

  
  // compute maximum buffer size
  k=0;
  for (i=0; i<BiD->nblocks; i++) {
      k=MAX(k,BiD->nblockcol[i]);
      k=MAX(k,BL->nblockrow[i]);
  } // end for i
  if (dbuff==NULL)
     buff=(FLOAT *)malloc((size_t)k*m*sizeof(FLOAT));
  else
     buff=dbuff;



  // --------------------------------------------------------------------------
  // forward substitution
  // pointer to the current diagonal block of rhssol
  psol=rhssol;
  // leading block diagonal index
  tt=0;
  for (i=0; i<BiD->nblocks; i++,psol+=k,tt+=k) {
      // size of the current diagonal block
      k=BiD->nblockcol[i];

      // downdate off-diagonal parts
      // auxiliary buffer
      pb=buff;
      // sub-diagonal index vector
      pi=BL->rowind[i];
      // number of sub-diagonal indices
      l =BL->nblockrow[i];
      // pointer to the sub-diagonal block
      pL=BL->valE[i];
      // scalar case, level-2-BLAS
      if (m==1 && l) {
	 // downdate sub-diagonal part
	 // scalar case, downdate in-place
	 if (k==1) {
	    // -diagonal entry
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    val=-*psol;
#else
	    val.r=-psol->r; val.i=-psol->i;
#endif // _SINGLE_REAL_ or _DOUBLE_REAL_

#ifdef _USE_MKL_
	    // use sparse AXPYI psol[pi]= val * pL + psol[pi]
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    CBLAS_AXPYI(l, val, pL,pi, psol-tt); 
#else
	    CBLAS_AXPYI(l, &val, pL,pi, psol-tt); 
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
	    
#else // hand-coded loop
	    // downshift
	    psol-=tt;
	    for (r=0; r<l; r++,pi++,pL++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	        psol[*pi]+=(*pL)*val;
#else
	        psol[*pi].r+=pL->r*val.r-pL->i*val.i;
	        psol[*pi].i+=pL->r*val.i+pL->i*val.r;
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_		
	    } // end for r
	    // upshift
	    psol+=tt;
#endif //-else _USE_MKL_
	    
	 }
	 else {
	    // gather entries from sol to buff
#ifdef _USE_MKL_
	    CBLAS_GTHR(l, psol-tt, buff,pi); 
#else // hand-coded loop
	    // downshift
	    psol-=tt;
	    for (r=0; r<l; r++,pi++)
	        *pb++=psol[*pi];
	    // upshift
	    psol+=tt;
#endif //-else _USE_MKL_
	    
	    // buff = -1*L21*sol + 1*buff
	    transa="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    alpha=-1.0;
	    beta = 1.0;
#else
	    alpha.r=-1.0; alpha.i=0.0;
	    beta.r = 1.0; beta.i =0.0;
#endif
	    i2=1;
	    GEMV(transa, &l,&k, &alpha, pL,&l,
		 psol,&i2, &beta,
		 buff,&i2, 1);
	    
	    // scatter entries from buff back to sol
#ifdef _USE_MKL_
	    CBLAS_SCTR(l, buff, pi, psol-tt); 
#else // hand-coded loop
	    // downshift
	    psol-=tt;
	    // sub-diagonal index vector
	    pi=BL->rowind[i];
	    // auxiliary buffer
	    pb=buff;
	    for (r=0; r<l; r++,pi++)
	        psol[*pi]=*pb++;
	    // upshift
	    psol+=tt;
#endif //-else _USE_MKL_
	    
	 } // end if-else k=1
      } // end if
      else if (l) { // level-3-BLAS case

	 // gather rows from sol to buff
	 // downshift
	 psol-=tt;
         for (r=0; r<l; r++,pi++,pb++)
	     // copy row *pi of sol to buff
	     COPY(&m, psol+*pi,&n, pb,&l);
	 // upshift
	 psol+=tt;
	 
	 // buff = -1*L21*sol + 1*buff
	 transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 alpha=-1.0;
	 beta = 1.0;
#else
	 alpha.r=-1.0; alpha.i=0.0;
	 beta.r = 1.0; beta.i =0.0;
#endif
	 i2=1;
	 GEMM(transa, transb, &l,&m,&k, &alpha, pL,&l,
              psol,&n, &beta,
              buff,&l, 1,1);
	 
         // scatter entries of buff back to sol
	 // downshift
	 psol-=tt;
	 // auxiliary buffer
	 pb=buff;
	 // sub-diagonal index vectos
	 pi=BL->rowind[i];
         for (r=0; r<l; r++,pi++,pb++)
	     // copy row *pi of buff to sol
	     COPY(&m, pb,&l, psol+*pi,&n);
	 // upshift
	 psol+=tt;
	 
      } // end if-else
  } // end for i
  // --------------------------------------------------------------------------

  
  
  // --------------------------------------------------------------------------
  // diagonal block multiplication
  // pointers to the current diagonal block of sol 
  psol=rhssol;
  // leading block diagonal index
  tt=0;
  for (i=0; i<BiD->nblocks; i++,psol+=k,tt+=k) {
      // size of the diagonal block
      k=BiD->nblockcol[i];
      // cache sol part ass. with the diagonal block (no alias!!!)
      if (invert_blocks)
	 for (j=0; j<m; j++)
	     memcpy(buff+k*j, psol+n*j, (size_t)k*sizeof(FLOAT));
      // pointer to the diagonal block
      pD=BiD->valD[i];
      // scalar case, level-2-BLAS
      if (m==1) {
	 // scalar case, downdate in-place
	 if (k==1) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    *psol*=*pD;
#else
	    val=*psol;
	    psol->r=val.r*pD->r-val.i*pD->i;
	    psol->i=val.r*pD->i+val.i*pD->r;
#endif
	 }
	 else {
	    if (invert_blocks) {
	       // sol = 1*D11*buff + 0*sol
	       transa="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	       alpha=1.0;
	       beta =0.0;
#else
	       alpha.r=1.0; alpha.i=0.0;
	       beta.r =0.0; beta.i =0.0;
#endif
	       i2=1;
	       GEMV(transa, &k,&k, &alpha, pD,&k,
		    buff,&i2, &beta,
		    psol,&i2, 1);
	    } // end if invert_blocks
	    else {
	      // D11*D11^H * sol <- sol
	      if (spd) {
		 POTRS("l", &k, &m, pD,&k, psol,&n, &ierr,1);
		 if (ierr<0) {
		    printf("LAPACK's POTRS: %ld-th argument had an illegal value\n", -ierr);
		    return;
		 }
	      }
	      else {
		 MYSYTRS("l", &k, &m, pD,&k, pivots+tt, psol,&n, &ierr,1);
		 if (ierr<0) {
		    printf("LAPACK's MYSYTRS: %ld-th argument had an illegal value\n", -ierr);
		    return;
		 }
	      }
	    } // end if-else invert_blocks
	 } // end if-else k=1
      } // end if
      else { // level-3-BLAS case
	 if (invert_blocks) {
	    // sol = 1*D11*buff + 0*sol
	    transa="n"; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    alpha=1.0;
	    beta =0.0;
#else
	    alpha.r=1.0; alpha.i=0.0;
	    beta.r =0.0; beta.i =0.0;
#endif
	    i2=1;
	    GEMM(transa, transb, &k,&m,&k, &alpha, pD,&k,
		 buff,&k, &beta,
		 psol,&n, 1,1);
	 } // end if invert_blocks
	 else {
	    // D11*D11^H * sol <- sol
	    if (spd) {
	       POTRS("l", &k, &m, pD,&k, psol,&n, &ierr,1);
	       if (ierr<0) {
		  printf("LAPACK's POTRS: %ld-th argument had an illegal value\n", -ierr);
		  return ;
	       }
	    }
	    else {
	       MYSYTRS("l", &k, &m, pD,&k, pivots+tt, psol,&n, &ierr,1);
	       if (ierr<0) {
		  printf("LAPACK's MYSYTRS: %ld-th argument had an illegal value\n", -ierr);
		  return ;
	       }
	    } // end if-else spd
	 } // end if-else invert_blocks
      } // end if-else
  } // end for i
  // --------------------------------------------------------------------------


  
  // --------------------------------------------------------------------------
  // backward substitution
  // pointers to the current diagonal block of sol
  psol=rhssol+n;
  // leading block diagonal index
  tt=n;
  for (i=BiD->nblocks-1; i>=0; i--) {
      // size of the diagonal block
      k=BiD->nblockcol[i];
      // adjust pointer and index w.r.t. the current block
      psol-=k; tt-=k;

      // downdate off-diagonal parts
      // auxiliary buffer
      pb=buff;
      // super-diagonal index vector
      pi=BL->rowind[i];
      // number of super-diagonal indices
      l =BL->nblockrow[i];
      // pointer to the transposed super-diagonal block
      pL=BL->valE[i];
      // scalar case, level-2-BLAS
      if (m==1 && l) {
	 // scalar case, downdate in-place
	 if (k==1) {

#ifdef _USE_MKL_
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    *psol-=CBLAS_DOTI(l, pL, pi, psol-tt);
#else
#ifdef _COMPLEX_SYMMETRIC_
	    CBLAS_DOTUI(l, pL, pi, psol-tt, &val);
#else
	    CBLAS_DOTI(l, pL, pi, psol-tt, &val);
#endif
	    psol->r-=val.r; psol->i-=val.i;
#endif //-else _SINGLE_REAL_ or _DOUBLE_REAL_
	   
#else // hand-coded loop   
	    // value to collect scalar product
	    val=*psol;
	    // downshift
	    psol-=tt;
	    // downdate by scalar product
	    for (r=0; r<l; r++,pi++,pL++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	        val-=(*pL)*psol[*pi];
#else
	        val.r-=pL->r*psol[*pi].r-CONJG(pL->i)*psol[*pi].i;
	        val.i-=pL->r*psol[*pi].i+CONJG(pL->i)*psol[*pi].r;
#endif
	    } // end for r
	    // upshift
	    psol+=tt;
	    *psol=val;
#endif //-else _USE_MKL_
	    
	 }
	 else {
	    // gather entries from sol to buff
#ifdef _USE_MKL_
	    CBLAS_GTHR(l, psol-tt, buff,pi); 
#else // hand-coded loop
	    // downshift
	    psol-=tt;
	    for (r=0; r<l; r++,pi++)
	        *pb++=psol[*pi];
	    // upshift
	    psol+=tt;
#endif //-else _USE_MKL_
	    
	    // sol = -1*BL21^T*buff + 1*sol
	    transa=TRANSA;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	    alpha=-1.0;
	    beta = 1.0;
#else
	    alpha.r=-1.0; alpha.i=0.0;
	    beta.r = 1.0; beta.i =0.0;
#endif
	    i2=1;
	    GEMV(transa, &l,&k, &alpha, pL,&l,
		 buff,&i2, &beta,
		 psol,&i2, 1);
	    
	 } // end if-else k=1
      } // end if
      else if (l) { // level-3-BLAS case
         // gather rows from sol to buff
	 // downshift
	 psol-=tt;
         for (r=0; r<l; r++,pi++,pb++)
	     // copy row *pi of sol to buff
	     COPY(&m, psol+*pi,&n, pb,&l);
	 // upshift
	 psol+=tt;
	 
	 // sol = -1*BL21^T*buff + 1*sol
	 transa=TRANSA; transb="n";
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 alpha=-1.0;
	 beta = 1.0;
#else
	 alpha.r=-1.0; alpha.i=0.0;
	 beta.r = 1.0; beta.i =0.0;
#endif
	 i2=1;
	 GEMM(transa, transb, &k,&m,&l, &alpha, pL,&l,
              buff,&l, &beta,
              psol,&n, 1,1);
	 
      } // end if-else
  } // end for i

  if (dbuff==NULL)
     free(buff);
} // end bildlsol
