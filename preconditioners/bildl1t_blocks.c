/* $Id: bildl1t_blocks.c 4078 2018-03-16 15:39:13Z bolle $

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

#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYBILDL1T_BLOCKS     SSYMbildl1t_blocks
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYBILDL1T_BLOCKS     DSYMbildl1t_blocks
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYBILDL1T_BLOCKS     CSYMbildl1t_blocks
#define CONJG(A)       (A)
#else // double complex
#define MYBILDL1T_BLOCKS     ZSYMbildl1t_blocks
#define CONJG(A)       (A)
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYBILDL1T_BLOCKS     SSYMbildl1t_blocks
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYBILDL1T_BLOCKS     DSYMbildl1t_blocks
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYBILDL1T_BLOCKS     CHERbildl1t_blocks
#define CONJG(A)       (-(A))
#else // double complex
#define MYBILDL1T_BLOCKS     ZHERbildl1t_blocks
#define CONJG(A)       (-(A))
#endif //-if-elif-else single-real

#endif //-else complex-symmetric



#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define ELBOW_BLOCK     1.3333
#define BLOCK_EXT       4
#define TWO_BY_TWO_THRESHOLD 0.1
#define TWO_BY_TWO_BOUND     1.5

// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

integer MYBILDL1T_BLOCKS(SPARSEMATRIX *A,
			 REALS *SL, REALS *SR, integer *p, integer *invq,
			 integer *blocksize, integer *nblocks,
			 REALS droptol)
{  
  /*
    given a reordering and scaling determine approximately a block structure
    by locally simulating ILDL(1,droptol), i.e., we only allow fill from the
    neighbouring original nodes and drop small entries in modulus


    Input
    -----
    A           assumed to be a nonsingular sparse matrix 
    p,invq      permutation p and inverse permutation invq associated with the a priori
                permutation matrix P
		we further assume that the row indices of A(p,p) in each column k=p(i)
		are already sorted in increasing order. 
    SL          real diagonal scaling matrix
    SR          not referenced

    Output
    ------
    blocksize   array with the initial block sizes of each diagonal block
                of the scaled and reordered system.
    nblocks     number of blocks within the array "blocksize"
  */

  integer n=A->nc,      // total size
          i,j,k,l,m,r,s,t, ii,jj,kk,tt,ss, // counters
          *pi,*pi2,    // index array pointers
    
          *buff,       // auxiliary index buffer
          nbuff,       // logical size of buff

          cnt_fill,    // counter for additional fill when merging blocks
          old_fill,    // fill without merging before/after considering step j
          add_fill,    // additional fill by step j without merging
          new_fill,    // total fill when merging with column/row j

          startblock,  // beginning and end of the current diagonal
          endblock,    // block
    
          *Afirst,     // physical position of the first nonzero entry of each 
                       // column as far as we still need to access this column
                       // to extract the associated rows
          *ATfirst,    // physical position of the first nonzero entry of each 
                       // row as far as we still need to access this row
                       // to extract the associated columns

    
          *idxL,       // temporary index list of row indices
          *idxposL,    // associated inverse mapping storing the position
          cntL,        // length of idxL
          cntL_old,    // old length of idxL when needed
          colj,colt,   // flags for indicating block column case
          flag,

          *idxLj,       // temporary index list of row indices column j
          *idxposLj,    // associated inverse mapping storing the position
          cntLj,        // length of idxLj
          cntLj_old;    // old length of idxLj when needed

  
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2,     // temporary numerical pointers
          *Adiag,       // array of diagonal entries 
          *Asdiag;      // array of sub-diagonal entries

  REALS   att,atj,ajt,ajtjt,aitit, // temporary matrix entries by modulus
          bnd,bndl,bndr,bndt, // bounds, parameters used for inverse of 2x2 blocks
          a11,a21,a22,   // auxiliary variables for determinant
          *pr,*pv,*pAT,  // temporary real numerical pointers
           rval;         // temporary scalar real numerical value
  DSparseMatrix AT;      // |SL*A(p,p)*SL|^T
  

  
  // create temporary index lists for rows and column
  idxL    =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposL =(integer *)calloc(n,sizeof(integer));
  // create additional temporary index lists at step j for rows and column
  idxLj   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposLj=(integer *)calloc(n,sizeof(integer));
  
  // buffer for index lists
  nbuff =n;
  buff  =(integer *)malloc((size_t)nbuff *sizeof(integer));

  // beginning of each column/row at step j
  Afirst =(integer *)calloc((size_t)n,sizeof(integer));
  ATfirst=(integer *)calloc((size_t)n,sizeof(integer));
  // diagonal entries 
  Adiag  =(FLOAT *)malloc(n*sizeof(FLOAT));
  // sub-diagonal entries 
  Asdiag =(FLOAT *)malloc(n*sizeof(FLOAT));
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
  for (i=0; i<n; i++)
      Adiag[i]=Asdiag[i]=0.0;
#else
  for (i=0; i<n; i++)
      Adiag[i].r=Asdiag[i].r=Adiag[i].i=Asdiag[i].i=0.0;
#endif

  // prepare memory for A^T, recall that for symmetry reasons, only half of A is
  // stored
  AT.nr=AT.nc=n;
  AT.ncol  =(integer *) calloc((size_t)n,sizeof(integer));
  AT.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
  AT.val   =(REALS   **)malloc((size_t)n*sizeof(REALS *));

  // count nnz for A^T
  pi=AT.ncol;
  for (i=0; i<n; i++) {
      // scan column A(p,p[i])
      // shortcuts
      ii=p[i]; 
      pi2=A->rowind[ii];
      l=A->ncol[ii];
      for (j=0; j<l; j++) {
	  // index k of A(p[k],p[i])
	  kk=pi2[j];
	  // p[k]=kk
	  k=invq[kk];
	  // increase number of elements in row k
	  (pi[k])++;
      } // end for j
  } // end for i

  // set up remaining memory for AT=CONJ(A(p,p)^T)
  for (i=0; i<n; i++) {
      // nz in row i of A(p,p)
      l=pi[i];
      AT.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
      AT.val[i]   =(REALS *)  malloc((size_t)l*sizeof(REALS));
      // reset AT.ncol
      pi[i]=0;
  } // end for i
#ifdef PRINT_INFO
  printf("BILDL1T: A^T set up (1.pass)\n");fflush(stdout);
#endif

  // copy: AT <- SL*CONJ(A(p,p)^T)*SL 
  // also scan A(p,p) in order to find the position of the diagonal entry
  // though it exists (otherwise use position that follows the diagonal entry)
  // and store it in Afirst
  if (SL==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 flag=-1;
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // A(p[k],p[i])
	     val=pA[j];
	     
	     // store the beginning of the diagonal/super-diagonal part
	     if (flag) {
	        if (k<i)
		   Afirst[i]=j+1;
		else {
		   Afirst[i]=j;
		   flag=0;
		} // end if-else
	     } // end if
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asdiag[i]=val;

	     // val <- CONJ(val)
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	     val.i=CONJG(val.i);
#endif
	     
	     // check for the super-diagonal entry A(p[i-1],p[i])=CONJ(A(p[i],p[i-1])
	     if (k==i-1)
	        Asdiag[i-1]=val;
	     
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     // AT(i,k) <- CONJ(A(p[k],p[i])
	     AT.val[k][m]   =FABS(val);
	     // increase nz in column k of AT
	     pi[k]=m+1;
	 } // end for j
     } // end for i
  }
  else { // SL!=0
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 flag=-1;
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // SL(p[k],p[k]) A(p[k],p[i]) SL(p[i],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     val=SL[kk]*pA[j]*SL[ii];
#else
	     val.r=SL[kk]*pA[j].r*SL[ii];
	     val.i=SL[kk]*pA[j].i*SL[ii];
#endif

	     // store the beginning of the diagonal/super-diagonal part
	     if (flag) {
	        if (k<i)
		   Afirst[i]=j+1;
		else {
		   Afirst[i]=j;
		   flag=0;
		} // end if-else
	     } // end if
	     // A(p[i],p[i]) found
	     if (k==i) 
	        Adiag[i]=val;
	     // check for existing sub-diagonal entry A(p[i+1],p[i])
	     else if (k==i+1)
	        Asdiag[i]=val;
	     
	     // val <- CONJ(val)
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
	     val.i=CONJG(val.i);
#endif
	     
	     // check for the super-diagonal entry A(p[i-1],p[i])=CONJ(A(p[i],p[i-1])
	     if (k==i-1)
	        Asdiag[i-1]=val;
	     
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     // AT(i,k) <- CONJ(A(p[k],p[i])
	     AT.val[k][m]   =FABS(val);
	     // increase nz in column k of AT
	     pi[k]=m+1;
	 } // end for j
     } // end for i
  } // end if-else SL=0
#ifdef PRINT_INFO
  printf("BILDL1T: A^T set up (2.pass)\n");fflush(stdout);
#endif


  // scan CONJ(A(p,p)^T) in order to find the position of the diagonal entry
  // though it exists (otherwise use position that follows the diagonal entry)
  // and store it in ATfirst
  for (i=0; i<n; i++) {
      // scan AT(:,i)=CONJ(A(p[i],p))
      // shortcuts
      pi2=AT.rowind[i];
      l=AT.ncol[i];
      for (j=0; j<l; j++) {
	  // index k of A(p[i],p[k])
	  k=pi2[j];
	  // diagonal and super-diagonal begins
	  if (k>=i)
	     break;
      } // end for j
      ATfirst[i]=j;
  } // end for i
#ifdef PRINT_INFO
  printf("BILDL1T: A^T set up (3.pass)\n");fflush(stdout);
#endif

#ifdef PRINT_INFO1
  printf("Afirst\n");
  for (i=0; i<n; i++)
      printf("%4ld",Afirst[i]);
  printf("\n");
  printf("ATfirst\n");
  for (i=0; i<n; i++)
      printf("%4ld",ATfirst[i]);
  printf("\n");
#endif
 
  // **************************************************************************
  // *****                       main loop                                *****
  // **************************************************************************
  cntL=0;   // counter for accumulated index lists
  cntLj=0;  // counter for temporary current index lists
  startblock=endblock=0; // initial beginning of the current block
  cnt_fill=0; // counter for additional fill-in when merging blocks
  *nblocks=0; // number of diagonal blocks of A
  for (j=0; j<n; j++) {

      // compute fill-in in column j of L
	    
#ifdef PRINT_INFO
      printf("startblock=%3ld, endblock=%3ld\n",startblock,endblock);
      for (m=0; m<cntL; m++)
	  printf("%6ld",idxL[m]);
      printf("\n");
      printf("column j=%3ld\n",j);
      fflush(stdout);
#endif

      // scalar bound
      bnd=FABS(Adiag[j]);
      // so far no additional column (i.e. no block approach)
      colj=-2;
      // check whether a 1x1 pivot or a 2x2 pivot should be chosen
      // first check preceding column j-1
      if (j>0) {
#ifdef PRINT_INFO
	 printf("also check column j-1=%3ld\n",j-1);
	 fflush(stdout);
#endif
	 // [a11 a21] = |A(p[j-1:j],p[j-1:j])|
	 // [a21 a22] 
	 a11=FABS(Adiag[j-1]);
	 a21=FABS(Asdiag[j-1]);
	 a22=bnd;
	 // check whether blocking the columns j-1:j is superior
	 // only take this into account, if the sub-diagonal entry is
	 // sufficiently large
 	 if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
	    // use 1/||A(p[j-1:j],p[j-1:j])^{-1}|| as measure
	    // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
	    val=Adiag[j-1]*Adiag[j]-Asdiag[j-1]*Asdiag[j-1];
#else
	    val.r=(Adiag[j-1].r *Adiag[j].r          -Adiag[j-1].i *Adiag[j].i)
	         -(Asdiag[j-1].r*Asdiag[j-1].r       -Asdiag[j-1].i*CONJG(Asdiag[j-1].i));
	    val.i=(Adiag[j-1].r *Adiag[j].i          +Adiag[j-1].i *Adiag[j].r)
	         -(Asdiag[j-1].r*CONJG(Asdiag[j-1].i)+Asdiag[j-1].i*Asdiag[j-1].r);
#endif
	    // 1/||A(p[j-1:j],p[j-1:j])^{-1}||
	    bndl=FABS(val)/MAX(a11+a21,a21+a22);

	    // check if a 2x2 pivot is preferred, otherwise skip it
	    if (bndl>TWO_BY_TWO_BOUND*a11) {
	       // use column j-1 as additional column
	       colj=j-1;
	       // printf("step %ld, block %ld:%ld, %12.4le,%12.4le\n",j,colj,j,bndl,a11);
	       // use as bound threshold/||A(p[j-1:j],p[j-1:j])^{-1}||
	       bnd=droptol*bndl;
	    } // end if
	    // else
	    //   printf("step %ld, bndl=%12.4le, bnd=%12.4le\n",j,bndl,a11);
	 } // end if a21 large enough
	 // else
	 //   printf("step %ld, A(p[%ld],p[%ld])=%12.4le is not large enough\n",j,j,j-1,a21);
      } // end if j>0
      
      // possibly second check the subsequent column j+1
      if (j<n-1 && colj<0) {
#ifdef PRINT_INFO
	 printf("also check column j+1=%3ld\n",j+1);
	 fflush(stdout);
#endif
	 // [a11 a21] = |A(p[j:j+1],p[j:j+1])|
	 // [a21 a22] 
	 a11=bnd;
	 a21=FABS(Asdiag[j]);
	 a22=FABS(Adiag[j+1]);
	 // check whether blocking columns j:j+1 is superior
	 // only take this into account, if the sub-diagonal entry is
	 // sufficiently large
 	 if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
	    // use 1/||A(p[j:j+1],p[j:j+1])^{-1}|| as measure
	    // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
	    val=Adiag[j]*Adiag[j+1]-Asdiag[j]*Asdiag[j];
#else
	    val.r=(Adiag[j].r *Adiag[j+1].r      -Adiag[j].i *Adiag[j+1].i)
	         -(Asdiag[j].r*Asdiag[j].r       -Asdiag[j].i*CONJG(Asdiag[j].i));
	    val.i=(Adiag[j].r *Adiag[j+1].i      +Adiag[j].i *Adiag[j+1].r)
	         -(Asdiag[j].r*CONJG(Asdiag[j].i)+Asdiag[j].i*Asdiag[j].r);
#endif
	    // 1/||A(p[j:j+1],p[j:j+1])^{-1}||
	    bndr=FABS(val)/MAX(a11+a21,a21+a22);

	    // check if a 2x2 pivot is preferred, otherwise skip it
	    if (bndr>TWO_BY_TWO_BOUND*a11) {
	       // use column j+1 as additional column
	       colj=j+1;
	       // printf("step %ld, block %ld:%ld, %12.4le,%12.4le\n",j,j,colj,bndr,a11);
	       // use as bound threshold/||A(p[j:j+1],p[j:j+1])^{-1}||
	       bnd=droptol*bndr;
	    } // end if
	    // else
	    //   printf("step %ld, bndr=%12.4le, bnd=%12.4le\n",j,bndr,a11);
	 } // end if a21 large enough
	 // else
	 //  printf("step %ld, A(p[%ld],p[%ld])=%12.4le is not large enough\n",j,j+1,j,a21);
      } // end if j<n-1

      // default bound is threshold*|A(p[j],p[j])|
      if (colj<0) {
	 // printf("step %ld, %12.4le\n",j,bnd);
	 bnd*=droptol;
      }
      // -----------------------------------------------
      // -----------------------------------------------
      // ------------ extract COLUMN j of A(p,p) -------
      // check sub-diagonal part of COLUMN j of A(p,p) for its numerical fill
#ifdef PRINT_INFO
      printf("extract column j=%3ld\n",j);
      fflush(stdout);
#endif
      // Afirst[j] starts with indices t>=j
      // shortcuts
      jj=p[j]; 
      pi=A->rowind[jj];
      pA=A->val[jj];
      l=A->ncol[jj];
      r=Afirst[j];
      for (m=r; m<l; m++) {
	  // row index tt of A(tt,p[j])
	  tt=pi[m];
	  // tt=p[t]
	  t=invq[tt];

	  // sub-diagonal fill
	  if ((colj<j&&t>j) || (colj>j&&t>j+1)) {
	     // extract associated numerical value
	     val=pA[m];
	     if (SL!=NULL) {
	        rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		val*=rval;
#else
		val.r*=rval;
		val.i*=rval;
#endif
	        rval=SL[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		val*=rval;
#else
		val.r*=rval;
		val.i*=rval;
#endif
	     } // end if

	     // check whether the entry is greater than a threshold relative to the diagonal entry
	     if (FABS(val)>=bnd) {
	        // a fill entry
	        idxLj[cntLj]=t;       // bookmark index and ...
		idxposLj[t]=++cntLj;  // its physical position, shifted by 1
	     } // end if 
	  } // end if t>j
      } // end for m
#ifdef PRINT_INFO1
      printf("L(:,%3ld): fill ILDL(0,%8.1le), only from A(p,p)\n",j,droptol);
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      // ------------ extract COLUMN j of AT=|SL*A(p,p)*SL|^T -----------
      // check sub-diagonal part of COLUMN j of |SL*A(p,p)*SL|^T for its
      // numerical fill
      // ATfirst[j] starts with indices t>=j
      // shortcuts
#ifdef PRINT_INFO
      printf("extract column j=%3ld of AT\n",j);
      fflush(stdout);
#endif
      pi=AT.rowind[j];
      pAT=AT.val[j];
      l=AT.ncol[j];
      r=ATfirst[j];
      // store counter before entering column j of AT
      cntLj_old=cntLj;
      for (m=r; m<l; m++) {
	  // row index t of AT(t,j)
	  t=pi[m];

	  // sub-diagonal fill
	  if ((colj<j&&t>j) || (colj>j&&t>j+1)) {
	     // check whether the entry is greater than a threshold relative to the diagonal entry
	     if (pAT[m]>=bnd) {
	        // a fill entry
	        idxLj[cntLj]=t;       // bookmark index and ...
		idxposLj[t]=++cntLj;  // its physical position, shifted by 1
	     } // end if 
	  } // end if t>j
      } // end for m

      // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
      // are distinct index arrays each with indices in increasing order
      // we temporarily merge them in buff if necessary
      if (0<cntLj_old && cntLj_old<cntLj) {
	 ii=0; jj=cntLj_old;
	 l=0;
	 while (ii<cntLj_old && jj<cntLj) {
	       t=idxLj[ii]; tt=idxLj[jj];
	       if (t<tt) {
		  buff[l++]=t;
		  ii++;
	       }
	       else { // tt<t, the case tt=t is avoided by flagging idxposLj
		  buff[l++]=tt;
		  jj++;
	       }
	 } // end while
	 while (ii<cntLj_old) {
	       t=idxLj[ii];
	       buff[l++]=t;
	       ii++;
	 } // end while
	 while (jj<cntLj) {
	       tt=idxLj[jj];
	       buff[l++]=tt;
	       jj++;
	 } // end while
	 // write the sorted list back to idxLj
	 for (l=0; l<cntLj; ) {
	     t=buff[l];
	     idxLj[l]=t;
	     idxposLj[t]=++l;
	 } // end for l
      } // end if 0<cntLj_old<cntLj


      // possibly we currently check a 2x2 pivot
      // additionally take column colj into account
      if (colj>=0) {
#ifdef PRINT_INFO
	 printf("check sub-diagonal part of column colj=%3ld of A\n",colj);
	 fflush(stdout);
#endif
	 // check sub-diagonal part of COLUMN colj of A(p,p) for its numerical fill
	 // Afirst[colj] starts with indices t>=j
	 // shortcuts
	 jj=p[colj]; 
	 pi=A->rowind[jj];
	 pA=A->val[jj];
	 l=A->ncol[jj];
	 r=Afirst[colj];
	 // store counter before entering column colj of A
	 cntLj_old=cntLj;
	 for (m=r; m<l; m++) {
	     // row index tt of A(tt,p[colj])
	     tt=pi[m];
	     // tt=p[t]
	     t=invq[tt];

	     // sub-diagonal fill, not previously considered by column j
	     if (!idxposLj[t]) {
	        if ((colj<j&&t>j) || (colj>j&&t>j+1)) {
		   // extract associated numerical value
		   val=pA[m];
		   if (SL!=NULL) {
		      rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		      val*=rval;
#else
		      val.r*=rval;
		      val.i*=rval;
#endif
		      rval=SL[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		      val*=rval;
#else
		      val.r*=rval;
		      val.i*=rval;
#endif
		   } // end if

		   // check whether the entry is greater than a threshold relative to the diagonal entry
		   if (FABS(val)>=bnd) {
		      // a fill entry
		      idxLj[cntLj]=t;       // bookmark index and ...
		      idxposLj[t]=++cntLj;  // its physical position, shifted by 1
		   } // end if 
		} // end if t>j
	     } // end if !idxposLj[t]
	 } // end for m
	 
	 // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	 // are distinct index arrays each with indices in increasing order
	 // we temporarily merge them in buff if necessary
	 if (0<cntLj_old && cntLj_old<cntLj) {
	    ii=0; jj=cntLj_old;
	    l=0;
	    while (ii<cntLj_old && jj<cntLj) {
	          t=idxLj[ii]; tt=idxLj[jj];
		  if (t<tt) {
		     buff[l++]=t;
		     ii++;
		  }
		  else { // tt<t, the case tt=t is avoided by flagging idxposLj
		     buff[l++]=tt;
		     jj++;
		  }
	    } // end while
	    while (ii<cntLj_old) {
	          t=idxLj[ii];
		  buff[l++]=t;
		  ii++;
	    } // end while
	    while (jj<cntLj) {
	          tt=idxLj[jj];
		  buff[l++]=tt;
		  jj++;
	    } // end while
	    // write the sorted list back to idxLj
	    for (l=0; l<cntLj; ) {
	        t=buff[l];
		idxLj[l]=t;
		idxposLj[t]=++l;
	    } // end for l
	 } // end if 0<cntLj_old<cntLj

	 // ------------ extract COLUMN colj of AT=|SL*A(p,p)*SL|^T ------------
	 // check sub-diagonal part of COLUMN colj of |SL*A(p,p)*SL|^T for its numerical fill
	 // ATfirst[colj] starts with indices t>=j
#ifdef PRINT_INFO
	 printf("check sub-diagonal part of column colj=%3ld of AT\n",colj);
	 fflush(stdout);
#endif
	 // shortcuts
	 pi=AT.rowind[colj];
	 pAT=AT.val[colj];
	 l=AT.ncol[colj];
	 r=ATfirst[colj];
	 // store counter before entering column colj of AT
	 cntLj_old=cntLj;
	 for (m=r; m<l; m++) {
	     // row index t of AT(t,colj)
	     t=pi[m];

	     // sub-diagonal fill, not previously considered by column j
	     if (!idxposLj[t]) {
	        if ((colj<j&&t>j) || (colj>j&&t>j+1)) {
		   // check whether the entry is greater than a threshold relative to the diagonal entry
		   if (pAT[m]>=bnd) {
		      // a fill entry
		      idxLj[cntLj]=t;       // bookmark index and ...
		      idxposLj[t]=++cntLj;  // its physical position, shifted by 1
		   } // end if 
		} // end if t>j
	     } // end if !idxposLj[t]
	 } // end for m
	 
	 // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	 // are distinct index arrays each with indices in increasing order
	 // we temporarily merge them in buff if necessary
	 if (0<cntLj_old && cntLj_old<cntLj) {
	    ii=0; jj=cntLj_old;
	    l=0;
	    while (ii<cntLj_old && jj<cntLj) {
	          t=idxLj[ii]; tt=idxLj[jj];
		  if (t<tt) {
		     buff[l++]=t;
		     ii++;
		  }
		  else { // tt<t, the case tt=t is avoided by flagging idxposLj
		     buff[l++]=tt;
		     jj++;
		  }
	    } // end while
	    while (ii<cntLj_old) {
	          t=idxLj[ii];
		  buff[l++]=t;
		  ii++;
	    } // end while
	    while (jj<cntLj) {
	          tt=idxLj[jj];
		  buff[l++]=tt;
		  jj++;
	    } // end while
	    // write the sorted list back to idxLj
	    for (l=0; l<cntLj; ) {
	        t=buff[l];
		idxLj[l]=t;
		idxposLj[t]=++l;
	    } // end for l
	 } // end if 0<cntLj_old<cntLj
      } // end if colj>=0
      
      
#ifdef PRINT_INFO
      if (colj<0)
	 printf("L(:,%3ld): fill ILDL(0,%8.1le) including AT=|SL*A(p,p)*SL|^T\n",
		j,droptol);
      else if (colj<j)
	 printf("L(:,%3ld:%3ld): fill ILDL(0,%8.1le) including AT=|SL*A(p,p)*SL|^T\n",
		colj,j,droptol);
      else // colj>j
	 printf("L(:,%3ld:%3ld): fill ILDL(0,%8.1le) including AT=|SL*A(p,p)*SL|^T\n",
		j,colj,droptol);
	
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      
      // now we know the entries of ILDL(0,droptol) w.r.t. L
      // next proceed with ILDL(1,droptol)
      // check COLUMNS t<j of A(p,p) for their numerical fill
      // the entries 0<=t<j are located before r=Afirst[j]
#ifdef PRINT_INFO
      printf("Afirst[%ld]=%ld\n",j,Afirst[j]);
#endif
      jj=p[j]; 
      pi=A->rowind[jj];
      pA=A->val[jj];
      r=Afirst[j];
      for (m=0; m<r && cntLj<(n-j-1); m++) {
	  // row index tt of A(tt,p[j])
	  tt=pi[m];
	  // tt=p[t]
	  t=invq[tt];
#ifdef PRINT_INFO
	  printf("t=%ld\n",t);
#endif

	  // super-diagonal entry, include column t of A(p,p) for ILDL(1,droptol)
	  if ((t<j&&colj!=j-1) || (t<j-1&&colj==j-1)) {
#ifdef PRINT_INFO
	     printf("... also check column t=%3ld\n",t);
	     fflush(stdout);
#endif
	     // check again whether to use a 2x2 pivot
	     // scalar bound
	     bndt=FABS(Adiag[t]);
	     // so far no additional column
	     colt=-2;
	     // check whether a 1x1 pivot or a 2x2 pivot should be chosen
	     // first check preceding column t-1
	     if (t>0) {
	        // [a11 a21] = |A(p[t-1:t],p[t-1:t])|
	        // [a21 a22] 
	        a11=FABS(Adiag[t-1]);
		a21=FABS(Asdiag[t-1]);
		a22=bndt;
	        // check whether blocking columns t-1:t is superior
	        // only take this into account, if the sub-diagonal entry is
	        // sufficiently large
		if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
		   // use 1/||A(p[t-1:t],p[t-1:t])^{-1}|| as measure
		   // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
		   val=Adiag[t-1]*Adiag[t]-Asdiag[t-1]*Asdiag[t-1];
#else
		   val.r=(Adiag[t-1].r *Adiag[t].r          -Adiag[t-1].i *Adiag[t].i)
		        -(Asdiag[t-1].r*Asdiag[t-1].r       -Asdiag[t-1].i*CONJG(Asdiag[t-1].i));
		   val.i=(Adiag[t-1].r *Adiag[t].i          +Adiag[t-1].i *Adiag[t].r)
		        -(Asdiag[t-1].r*CONJG(Asdiag[t-1].i)+Asdiag[t-1].i*Asdiag[t-1].r);
#endif
		   // 1/||A(p[t-1:t],p[t-1:t])^{-1}||
		   bndl=FABS(val)/MAX(a11+a21,a21+a22);

		   // check if a 2x2 pivot is preferred, otherwise skip it
		   if (bndl>TWO_BY_TWO_BOUND*a11) {
		      // use column t-1 as additional column
		      colt=t-1;
		      // use as bound 1/||A(p[t-1:t],p[t-1:t])^{-1}||
		      bndt=bndl;
		   } // end if
		} // end if a21 large enough
	     } // end if t>0
      
	     // second check the subsequent column t+1<j (resp. t<j-1)
	     if (colt<0 && ((colj==j-1&&t+1<colj) || (colj!=j-1&&t+1<j))) {
#ifdef PRINT_INFO
	        printf("... also check column t+1=%3ld\n",t+1);
		fflush(stdout);
#endif
	        // [a11 a21] = |A(p[t:t+1],p[t:t+1])|
	        // [a21 a22] 
	        a11=bndt;
		a21=FABS(Asdiag[t]);
		a22=FABS(Adiag[t+1]);
	        // check whether blocking columns t:t+1 is superior
	        // only take this into account, if the sub-diagonal entry is
	        // sufficiently large
		if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
		   // use 1/||A(p[t:t+1],p[t:t+1])^{-1}|| as measure
		   // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
		   val=Adiag[t]*Adiag[t+1]-Asdiag[t]*Asdiag[t];
#else
		   val.r=(Adiag[t].r *Adiag[t+1].r      -Adiag[t].i *Adiag[t+1].i)
		        -(Asdiag[t].r*Asdiag[t].r       -Asdiag[t].i*CONJG(Asdiag[t].i));
		   val.i=(Adiag[t].r *Adiag[t+1].i      +Adiag[t].i *Adiag[t+1].r)
		        -(Asdiag[t].r*CONJG(Asdiag[t].i)+Asdiag[t].i*Asdiag[t].r);
#endif
		   // 1/||A(p[t:t+1],p[t:t+1])^{-1}||
		   bndr=FABS(val)/MAX(a11+a21,a21+a22);

		   // check if a 2x2 pivot is preferred, otherwise skip it
		   if (bndr>TWO_BY_TWO_BOUND*a11) {
		      // use column t+1 as additional column
		      colt=t+1;
		      // use as bound 1/||A(p[t:t+1],p[t:t+1])^{-1}||
		      bndt=bndr;
		   } // end if
		} // end if a21 large enough
	     } // end if t+1<j ...

	     // default bound is bndt=|A(p[t],p[t])|
	     
	     // extract associated numerical value A(p[t],p[j])
	     val=pA[m];
	     if (SL!=NULL) {
	        rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		val*=rval;
#else
		val.r=rval;
		val.i=rval;
#endif
	        rval=SL[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
	        val*=rval;
#else
		val.r*=rval;
		val.i*=rval;
#endif
	     } // end if SL!=0
	     // |A(p[t],p[j])|
	     atj=FABS(val);
	     
	     //                                  threshold
	     // ---------------------------------------------------------------------------
	     // ||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||*||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
	     att=bndt*bnd;

	     // search for beginning of A(p[j:end],p[t])
	     // right now Afirst[t] is only guaranteed to start with some
	     // index s>=t
	     // shortcuts
	     tt=p[t]; 
	     pi2=A->rowind[tt];
	     pA2=A->val[tt];
	     ii=A->ncol[tt];
	     for (l=Afirst[t]; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         ss=pi2[l];
		 // ss=p[s]
		 s=invq[ss];
		 // trailing part
	         if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		    break;
	     } // end for l
	     Afirst[t]=l;
#ifdef PRINT_INFO
	     printf("Afirst[%ld]=%ld\n",t,Afirst[t]);
#endif
	     
	     // check COLUMN t<j of A(p,p) for fill-in
	     // right now, all t>j are located at position Afirst[t] and higher
	     // store counter before entering column t<j
	     cntLj_old=cntLj;
	     for (; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         ss=pi2[l];
		 // ss=p[s]
		 s=invq[ss];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposLj[s]) {	     
		    // extract associated numerical value A(p[s],p[t])
		    val=pA2[l];
		    if (SL!=NULL) {
		       rval=SL[ss];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		       rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		    } // end if SL!=0

		    // threshold relative to the inv. diagonal blocks
		    // |A(p[s],p[t])A(p[t],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		    if (FABS(val)*atj>=att) {
		       // a fill entry
		       idxLj[cntLj]=s;       // bookmark index and ...
		       idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l

	     // fill of column t of A(p,p) now incorporated

	     // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	     // are distinct index arrays each with indices in increasing order
	     // we temporarily merge them in buff if necessary
	     if (0<cntLj_old && cntLj_old<cntLj) {
	        ii=0; jj=cntLj_old;
		l=0;
		while (ii<cntLj_old && jj<cntLj) {
	              t=idxLj[ii]; tt=idxLj[jj];
		      if (t<tt) {
			 buff[l++]=t;
			 ii++;
		      }
		      else { // tt<t, the case tt=t is avoided by flagging idxposLj
			 buff[l++]=tt;
			 jj++;
		      }
		} // end while
		while (ii<cntLj_old) {
		      t=idxLj[ii];
		      buff[l++]=t;
		      ii++;
		} // end while
		while (jj<cntLj) {
		      tt=idxLj[jj];
		      buff[l++]=tt;
		      jj++;
		} // end while
		// write the sorted list back to idxLj
		for (l=0; l<cntLj; ) {
		    t=buff[l];
		    idxLj[l]=t;
		    idxposLj[t]=++l;
		} // end for l
	     } // end if 0<cntLj_old<cntLj


	     // search for beginning of A^T(p[j:end],p[t])
	     // right now ATfirst[t] is only guaranteed to start with some
	     // index s>=t
	     // shortcuts
	     pi2=AT.rowind[t];
	     pv =AT.val[t];
	     ii=AT.ncol[t];
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A^T(p[s],p[t])
	         s=pi2[l];
		 // trailing part
	         if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		    break;
	     } // end for l
	     ATfirst[t]=l;
#ifdef PRINT_INFO
	     printf("ATfirst[%ld]=%ld\n",t,ATfirst[t]);
#endif
	     
	     // check COLUMN t<j of A^T(p,p) for fill-in
	     // right now, all t>j are located at position ATfirst[t] and higher
	     // store counter before entering column t<j of AT
	     cntLj_old=cntLj;
	     for (; l<ii; l++) {
	         // row index s of A^T(p[s],p[t])
	         s=pi2[l];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposLj[s]) {	     
		    // threshold relative to the inv. diagonal blocks
		    // |A^T(p[s],p[t])A(p[t],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		    if (pv[l]*atj>=att) {
		       // a fill entry
		       idxLj[cntLj]=s;       // bookmark index and ...
		       idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l
	     
	     // fill of column t of AT now incorporated
	     
	     // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	     // are distinct index arrays each with indices in increasing order
	     // we temporarily merge them in buff if necessary
	     if (0<cntLj_old && cntLj_old<cntLj) {
	        ii=0; jj=cntLj_old;
		l=0;
		while (ii<cntLj_old && jj<cntLj) {
	              t=idxLj[ii]; tt=idxLj[jj];
		      if (t<tt) {
			 buff[l++]=t;
			 ii++;
		      }
		      else { // tt<t, the case tt=t is avoided by flagging idxposLj
			 buff[l++]=tt;
			 jj++;
		      }
		} // end while
		while (ii<cntLj_old) {
		      t=idxLj[ii];
		      buff[l++]=t;
		      ii++;
		} // end while
		while (jj<cntLj) {
		      tt=idxLj[jj];
		      buff[l++]=tt;
		      jj++;
		} // end while
		// write the sorted list back to idxLj
		for (l=0; l<cntLj; ) {
		    t=buff[l];
		    idxLj[l]=t;
		    idxposLj[t]=++l;
		} // end for l
	     } // end if 0<cntLj_old<cntLj

	     // possibly consider fill by column colt
	     if (colt>=0) {
#ifdef PRINT_INFO
	        printf("possibly consider fill by column colt=%3ld\n",colt);
		fflush(stdout);
#endif
	        // search for beginning of A(p[j:end],p[colt])
	        // right now Afirst[t] is only guaranteed to start with some
	        // index s>=colt
	        // shortcuts
	        tt=p[colt]; 
		pi2=A->rowind[tt];
		pA2=A->val[tt];
		ii=A->ncol[tt];
		for (l=Afirst[colt]; l<ii; l++) {
	            // row index ss of A(ss,p[colt])
	            ss=pi2[l];
		    // ss=p[s]
		    s=invq[ss];
		    // trailing part
		    if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		       break;
		} // end for l
		Afirst[colt]=l;
#ifdef PRINT_INFO
		printf("Afirst[%ld]=%ld\n",colt,Afirst[colt]);
#endif
	     
		// check COLUMN colt<j of A(p,p) for fill-in
		// right now, all colt>j are located at position Afirst[colt] and higher
		// store counter before entering column colt<j
		cntLj_old=cntLj;
		for (; l<ii; l++) {
		    // row index ss of A(ss,p[colt])
	            ss=pi2[l];
		    // ss=p[s]
		    s=invq[ss];
	  
		    // sub-diagonal entry that has not been considered as fill yet
		    if (!idxposLj[s]) {	     
		       // extract associated numerical value A(p[s],p[colt])
		       val=pA2[l];
		       if (SL!=NULL) {
		          rval=SL[ss];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			  val*=rval;
#else
			  val.r*=rval;
			  val.i*=rval;
#endif
			  rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			  val*=rval;
#else
			  val.r*=rval;
			  val.i*=rval;
#endif
		       } // end if SL!=0

		       // threshold relative to the diagonal entry
		       // |A(p[s],p[colt])A(p[colt],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		       if (FABS(val)*atj>=att) {
			  // a fill entry
			  idxLj[cntLj]=s;       // bookmark index and ...
			  idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		       } // end if 
		    } // end if s>j
		} // end for l

		// fill of column colt of A(p,p) now incorporated

		// now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
		// are distinct index arrays each with indices in increasing order
		// we temporarily merge them in buff if necessary
		if (0<cntLj_old && cntLj_old<cntLj) {
		   ii=0; jj=cntLj_old;
		   l=0;
		   while (ii<cntLj_old && jj<cntLj) {
	                 t=idxLj[ii]; tt=idxLj[jj];
			 if (t<tt) {
			    buff[l++]=t;
			    ii++;
			 }
			 else { // tt<t, the case tt=t is avoided by flagging idxposLj
			    buff[l++]=tt;
			    jj++;
			 }
		   } // end while
		   while (ii<cntLj_old) {
		         t=idxLj[ii];
			 buff[l++]=t;
			 ii++;
		   } // end while
		   while (jj<cntLj) {
		         tt=idxLj[jj];
			 buff[l++]=tt;
			 jj++;
		   } // end while
		   // write the sorted list back to idxLj
		   for (l=0; l<cntLj; ) {
		       t=buff[l];
		       idxLj[l]=t;
		       idxposLj[t]=++l;
		   } // end for l
		} // end if 0<cntLj_old<cntLj


		// search for beginning of A^T(p[j:end],p[colt])
		// right now ATfirst[colt] is only guaranteed to start with some
		// index s>=colt
		// shortcuts
		pi2=AT.rowind[colt];
		pv =AT.val[colt];
		ii=AT.ncol[colt];
		for (l=ATfirst[colt]; l<ii; l++) {
	            // row index ss of A^T(p[s],p[colt])
		    s=pi2[l];
		    // trailing part
		    if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		       break;
		} // end for l
		ATfirst[colt]=l;
#ifdef PRINT_INFO
		printf("ATfirst[%ld]=%ld\n",colt,ATfirst[colt]);
#endif
	     
		// check COLUMN colt<j of A^T(p,p) for fill-in
		// right now, all colt>j are located at position ATfirst[colt] and higher
		// store counter before entering column colt<j of AT
		cntLj_old=cntLj;
		for (; l<ii; l++) {
	            // row index s of A^T(p[s],p[colt])
	            s=pi2[l];
	  
		    // sub-diagonal entry that has not been considered as fill yet
		    if (!idxposLj[s]) {	     
		       // threshold relative to the inv. diagonal block
		       // |A^T(p[s],p[colt])A(p[colt],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		       if (pv[l]*atj>=att) {
		 	  // a fill entry
			  idxLj[cntLj]=s;       // bookmark index and ...
			  idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		       } // end if 
		    } // end if s>j
		} // end for l
	     
		// fill of column colt of AT now incorporated
	     
		// now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
		// are distinct index arrays each with indices in increasing order
		// we temporarily merge them in buff if necessary
		if (0<cntLj_old && cntLj_old<cntLj) {
		   ii=0; jj=cntLj_old;
		   l=0;
		   while (ii<cntLj_old && jj<cntLj) {
	                 t=idxLj[ii]; tt=idxLj[jj];
			 if (t<tt) {
			    buff[l++]=t;
			    ii++;
			 }
			 else { // tt<t, the case tt=t is avoided by flagging idxposLj
			    buff[l++]=tt;
			    jj++;
			 }
		   } // end while
		   while (ii<cntLj_old) {
		         t=idxLj[ii];
			 buff[l++]=t;
			 ii++;
		   } // end while
		   while (jj<cntLj) {
		         tt=idxLj[jj];
			 buff[l++]=tt;
			 jj++;
		   } // end while
		   // write the sorted list back to idxLj
		   for (l=0; l<cntLj; ) {
		       t=buff[l];
		       idxLj[l]=t;
		       idxposLj[t]=++l;
		   } // end for l
		} // end if 0<cntLj_old<cntLj
	     } // end if colt>=0
	  } // end if t<j
      } // end for m

      // check COLUMNS t<j of A^T(p,p) for their numerical fill
      // the entries 0<=t<j are located before r=ATfirst[j]
#ifdef PRINT_INFO
      printf("ATfirst[%ld]=%ld\n",j,ATfirst[j]);
#endif
      pi=AT.rowind[j];
      pAT=AT.val[j];
      r=ATfirst[j];
      for (m=0; m<r && cntLj<(n-j-1); m++) {
	  // row index t of A(p[t],p[j])
	  t=pi[m];
#ifdef PRINT_INFO
	  printf("t=%ld\n",t);
#endif

	  // super-diagonal entry, include column t of A^T(p,p) for ILDL(1,droptol)
	  if ((t<j&&colj!=j-1) || (t<j-1&&colj==j-1)) {
#ifdef PRINT_INFO
	     printf("... also check column t=%3ld\n",t);
	     fflush(stdout);
#endif
	     // check again whether to use a 2x2 pivot
	     // scalar bound
	     bndt=FABS(Adiag[t]);
	     // so far no additional column
	     colt=-2;
	     // check whether a 1x1 pivot or a 2x2 pivot should be chosen
	     // first check preceding column t-1
	     if (t>0) {
	        // [a11 a21] = |A(p[t-1:t],p[t-1:t])|
	        // [a21 a22] 
	        a11=FABS(Adiag[t-1]);
		a21=FABS(Asdiag[t-1]);
		a22=bndt;
	        // check whether blocking column t-1:t is superior
	        // only take this into account, if the sub-diagonal entry is
	        // sufficiently large
		if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
		   // use 1/||A(p[t-1:t],p[t-1:t])^{-1}|| as measure
		   // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
		   val=Adiag[t-1]*Adiag[t]-Asdiag[t-1]*Asdiag[t-1];
#else
		   val.r=(Adiag[t-1].r *Adiag[t].r          -Adiag[t-1].i *Adiag[t].i)
		        -(Asdiag[t-1].r*Asdiag[t-1].r       -Asdiag[t-1].i*CONJG(Asdiag[t-1].i));
		   val.i=(Adiag[t-1].r *Adiag[t].i          +Adiag[t-1].i *Adiag[t].r)
		        -(Asdiag[t-1].r*CONJG(Asdiag[t-1].i)+Asdiag[t-1].i*Asdiag[t-1].r);
#endif
		   // 1/||A(p[t-1:t],p[t-1:t])^{-1}||
		   bndl=FABS(val)/MAX(a11+a21,a21+a22);

		   // check if a 2x2 pivot is preferred, otherwise skip it
		   if (bndl>TWO_BY_TWO_BOUND*a11) {
		      // use column t-1 as additional column
		      colt=t-1;
		      // use as bound 1/||A(p[t-1:t],p[t-1:t])^{-1}||
		      bndt=bndl;
		   } // end if
		} // end if a21 large enough
	     } // end if t>0
      
	     // second check the subsequent column t+1<j (resp. t<j-1)
	     if (colt<0 && ((colj==j-1&&t+1<colj) || (colj!=j-1&&t+1<j))) {
	        // [a11 a21] = |A(p[t:t+1],p[t:t+1])|
	        // [a21 a22] 
	        a11=bndt;
		a21=FABS(Asdiag[t]);
		a22=FABS(Adiag[t+1]);
	        // check whether blocking columns t:t+1 is superior
	        // only take this into account, if the sub-diagonal entry is
	        // sufficiently large
		if (a21>TWO_BY_TWO_THRESHOLD*MIN(a11,a22)) {
		   // use 1/||A(p[t:t+1],p[t:t+1])^{-1}|| as measure
		   // 2x2 determinant, modulus
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_	
		   val=Adiag[t]*Adiag[t+1]-Asdiag[t]*Asdiag[t];
#else
		   val.r=(Adiag[t].r *Adiag[t+1].r      -Adiag[t].i *Adiag[t+1].i)
		        -(Asdiag[t].r*Asdiag[t].r       -Asdiag[t].i*CONJG(Asdiag[t].i));
		   val.i=(Adiag[t].r *Adiag[t+1].i      +Adiag[t].i *Adiag[t+1].r)
		        -(Asdiag[t].r*CONJG(Asdiag[t].i)+Asdiag[t].i*Asdiag[t].r);
#endif
		   // 1/||A(p[t:t+1],p[t:t+1])^{-1}||
		   bndr=FABS(val)/MAX(a11+a21,a21+a22);

		   // check if a 2x2 pivot is preferred, otherwise skip it
		   if (bndr>TWO_BY_TWO_BOUND*a11) {
		      // use column t+1 as additional column
		      colt=t+1;
		      // use as bound 1/||A(p[t:t+1],p[t:t+1])^{-1}||
		      bndt=bndr;
		   } // end if
		} // end if a21 large enough
	     } // end if t+1<j ...

	     // default bound is |A(p[t],p[t])|
	     
	     // extract associated numerical value A^T(p[t],p[j])
	     // |A^T(p[t],p[j])|
	     atj=pAT[m];
	     
	     //                                  threshold
	     // ---------------------------------------------------------------------------
	     // ||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||*||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
	     att=bndt*bnd;

	     // search for beginning of A(p[j:end],p[t])
	     // right now Afirst[t] is only guaranteed to start with some
	     // index s>=t
	     // shortcuts
	     tt=p[t]; 
	     pi2=A->rowind[tt];
	     pA2=A->val[tt];
	     ii=A->ncol[tt];
	     for (l=Afirst[t]; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         ss=pi2[l];
		 // ss=p[s]
		 s=invq[ss];
		 // trailing part
	         if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		    break;
	     } // end for l
	     Afirst[t]=l;
#ifdef PRINT_INFO
	     printf("Afirst[%ld]=%ld\n",t,Afirst[t]);
#endif
	     
	     // check COLUMN t<j of A(p,p) for fill-in
	     // right now, all t>j are located at position Afirst[t] and higher
	     // store counter before entering column t<j
	     cntLj_old=cntLj;
	     for (; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         ss=pi2[l];
		 // ss=p[s]
		 s=invq[ss];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposLj[s]) {	     
		    // extract associated numerical value A(p[s],p[t])
		    val=pA2[l];
		    if (SL!=NULL) {
		       rval=SL[ss];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		       rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		    } // end if SL!=0

		    // threshold relative to the inv. diagonal blocks
		    // |A(p[s],p[t])A(p[t],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		    if (FABS(val)*atj>=att) {
		       // a fill entry
		       idxLj[cntLj]=s;       // bookmark index and ...
		       idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l

	     // fill of column t of A(p,p) now incorporated

	     // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	     // are distinct index arrays each with indices in increasing order
	     // we temporarily merge them in buff if necessary
	     if (0<cntLj_old && cntLj_old<cntLj) {
	        ii=0; jj=cntLj_old;
		l=0;
		while (ii<cntLj_old && jj<cntLj) {
	              t=idxLj[ii]; tt=idxLj[jj];
		      if (t<tt) {
			 buff[l++]=t;
			 ii++;
		      }
		      else { // tt<t, the case tt=t is avoided by flagging idxposLj
			 buff[l++]=tt;
			 jj++;
		      }
		} // end while
		while (ii<cntLj_old) {
		      t=idxLj[ii];
		      buff[l++]=t;
		      ii++;
		} // end while
		while (jj<cntLj) {
		      tt=idxLj[jj];
		      buff[l++]=tt;
		      jj++;
		} // end while
		// write the sorted list back to idxLj
		for (l=0; l<cntLj; ) {
		    t=buff[l];
		    idxLj[l]=t;
		    idxposLj[t]=++l;
		} // end for l
	     } // end if 0<cntLj_old<cntLj


	     // search for beginning of A^T(p[j:end],p[t])
	     // right now ATfirst[t] is only guaranteed to start with some
	     // index s>=t
	     // shortcuts
	     pi2=AT.rowind[t];
	     pv=AT.val[t];
	     ii=AT.ncol[t];
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A^T(p[s],p[t])
	         s=pi2[l];
		 // trailing part
	         if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		    break;
	     } // end for l
	     ATfirst[t]=l;
#ifdef PRINT_INFO
	     printf("ATfirst[%ld]=%ld\n",t,ATfirst[t]);
#endif
	     
	     // check COLUMN t<j of A^T(p,p) for fill-in
	     // right now, all t>j are located at position Afirst[t] and higher
	     // store counter before entering column t<j of AT
	     cntLj_old=cntLj;
	     for (; l<ii; l++) {
	         // row index s of A(p[s],p[t])
	         s=pi2[l];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposLj[s]) {	     
		    // threshold relative to the inv. diagonal block
		    // |A^T(p[s],p[t])A(p[t],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		    if (pv[l]*atj>=att) {
		       // a fill entry
		       idxLj[cntLj]=s;       // bookmark index and ...
		       idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l
	     
	     // fill of column t of AT now incorporated
	     
	     // now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
	     // are distinct index arrays each with indices in increasing order
	     // we temporarily merge them in buff if necessary
	     if (0<cntLj_old && cntLj_old<cntLj) {
	        ii=0; jj=cntLj_old;
		l=0;
		while (ii<cntLj_old && jj<cntLj) {
	              t=idxLj[ii]; tt=idxLj[jj];
		      if (t<tt) {
			 buff[l++]=t;
			 ii++;
		      }
		      else { // tt<t, the case tt=t is avoided by flagging idxposLj
			 buff[l++]=tt;
			 jj++;
		      }
		} // end while
		while (ii<cntLj_old) {
		      t=idxLj[ii];
		      buff[l++]=t;
		      ii++;
		} // end while
		while (jj<cntLj) {
		      tt=idxLj[jj];
		      buff[l++]=tt;
		      jj++;
		} // end while
		// write the sorted list back to idxLj
		for (l=0; l<cntLj; ) {
		    t=buff[l];
		    idxLj[l]=t;
		    idxposLj[t]=++l;
		} // end for l
	     } // end if 0<cntLj_old<cntLj

	     // possibly consider fill by column colt
	     if (colt>=0) {
	        // search for beginning of A(p[j:end],p[colt])
	        // right now Afirst[colt] is only guaranteed to start with some
	        // index s>=colt
	        // shortcuts
	        tt=p[colt]; 
		pi2=A->rowind[tt];
		pA2=A->val[tt];
		ii=A->ncol[tt];
		for (l=Afirst[colt]; l<ii; l++) {
		    // row index ss of A(ss,p[t])
		    ss=pi2[l];
		    // ss=p[s]
		    s=invq[ss];
		    // trailing part
		    if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		       break;
		} // end for l
		Afirst[colt]=l;
#ifdef PRINT_INFO
		printf("Afirst[%ld]=%ld\n",colt,Afirst[colt]);
#endif
	     
		// check COLUMN colt<j of A(p,p) for fill-in
		// right now, all colt>j are located at position Afirst[colt] and higher
		// store counter before entering column colt<j
		cntLj_old=cntLj;
		for (; l<ii; l++) {
	            // row index ss of A(ss,p[t])
		    ss=pi2[l];
		    // ss=p[s]
		    s=invq[ss];
	  
		    // sub-diagonal entry that has not been considered as fill yet
		    if (!idxposLj[s]) {	     
		       // extract associated numerical value A(p[s],p[t])
		       val=pA2[l];
		       if (SL!=NULL) {
			  rval=SL[ss];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			  val*=rval;
#else
			  val.r*=rval;
			  val.i*=rval;
#endif
			  rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
			  val*=rval;
#else
			  val.r*=rval;
			  val.i*=rval;
#endif
		       } // end if SL!=0

		       // threshold relative to the inv. diagonal blocks
		       // |A(p[s],p[colt])A(p[colt],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		       if (FABS(val)*atj>=att) {
			  // a fill entry
			  idxLj[cntLj]=s;       // bookmark index and ...
			  idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		       } // end if 
		    } // end if s>j
		} // end for l

		// fill of column colt of A(p,p) now incorporated

		// now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
		// are distinct index arrays each with indices in increasing order
		// we temporarily merge them in buff if necessary
		if (0<cntLj_old && cntLj_old<cntLj) {
		   ii=0; jj=cntLj_old;
		   l=0;
		   while (ii<cntLj_old && jj<cntLj) {
	                 t=idxLj[ii]; tt=idxLj[jj];
			 if (t<tt) {
			    buff[l++]=t;
			    ii++;
			 }
			 else { // tt<t, the case tt=t is avoided by flagging idxposLj
			    buff[l++]=tt;
			    jj++;
			 }
		   } // end while
		   while (ii<cntLj_old) {
		         t=idxLj[ii];
			 buff[l++]=t;
			 ii++;
		   } // end while
		   while (jj<cntLj) {
		         tt=idxLj[jj];
			 buff[l++]=tt;
			 jj++;
		   } // end while
		   // write the sorted list back to idxLj
		   for (l=0; l<cntLj; ) {
		       t=buff[l];
		       idxLj[l]=t;
		       idxposLj[t]=++l;
		   } // end for l
		} // end if 0<cntLj_old<cntLj


		// search for beginning of A^T(p[j:end],p[colt])
		// right now ATfirst[colt] is only guaranteed to start with some
		// index s>=colt
		// shortcuts
		pi2=AT.rowind[colt];
		pv=AT.val[colt];
		ii=AT.ncol[colt];
		for (l=ATfirst[colt]; l<ii; l++) {
		    // row index ss of A^T(p[s],p[colt])
		    s=pi2[l];
		    // trailing part
		    if ((colj<j&&s>j) || (colj>j&&s>j+1)) 
		       break;
		} // end for l
		ATfirst[colt]=l;
#ifdef PRINT_INFO
		printf("ATfirst[%ld]=%ld\n",colt,ATfirst[colt]);
#endif
	     
		// check COLUMN colt<j of A^T(p,p) for fill-in
		// right now, all colt>j are located at position Afirst[colt] and higher
		// store counter before entering column colt<j of AT
		cntLj_old=cntLj;
		for (; l<ii; l++) {
		    // row index s of A(p[s],p[colt])
		    s=pi2[l];
	  
		    // sub-diagonal entry that has not been considered as fill yet
		    if (!idxposLj[s]) {	     
		       // threshold relative to the inv. diagonal block
		       // |A^T(p[s],p[colt])A(p[colt],p[j])| >= droptol/||A(p[t-1,t,t+1],p[t-1,t,t+1])^{-1}||/||A(p[j-1,j,j+1],p[j-1,j,j+1])^{-1}||
		       if (pv[l]*atj>=att) {
			  // a fill entry
			  idxLj[cntLj]=s;       // bookmark index and ...
			  idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		       } // end if 
		    } // end if s>j
		} // end for l
	     
		// fill of column colt of AT now incorporated
	     
		// now idxLj[0,...,cntLj_old-1] and idxLj[cntLj_old,...,cntLj-1]
		// are distinct index arrays each with indices in increasing order
		// we temporarily merge them in buff if necessary
		if (0<cntLj_old && cntLj_old<cntLj) {
		   ii=0; jj=cntLj_old;
		   l=0;
		   while (ii<cntLj_old && jj<cntLj) {
		         t=idxLj[ii]; tt=idxLj[jj];
			 if (t<tt) {
			    buff[l++]=t;
			    ii++;
			 }
			 else { // tt<t, the case tt=t is avoided by flagging idxposLj
			    buff[l++]=tt;
			    jj++;
			 }
		   } // end while
		   while (ii<cntLj_old) {
		         t=idxLj[ii];
			 buff[l++]=t;
			 ii++;
		   } // end while
		   while (jj<cntLj) {
		         tt=idxLj[jj];
			 buff[l++]=tt;
			 jj++;
		   } // end while
		   // write the sorted list back to idxLj
		   for (l=0; l<cntLj; ) {
		       t=buff[l];
		       idxLj[l]=t;
		       idxposLj[t]=++l;
		   } // end for l
		} // end if 0<cntLj_old<cntLj
	     } // end if colt>=0
	  } // end if t<j
      } // end for m

      
#ifdef PRINT_INFO
      if (colj<0)
	 printf("L(:,%3ld): fill ILDL(1,%8.1le)\n",j,droptol);
      else if (colj==j-1)
	 printf("L(:,%3ld:%3ld): fill ILDL(1,%8.1le)\n",j-1,j,droptol);
      else // colj==j+1
	 printf("L(:,%3ld:%3ld): fill ILDL(1,%8.1le)\n",j,j+1,droptol);
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      // -----------------------------------------------
      // -----------------------------------------------
      // ---------- END extract COLUMN j of A ----------

      
      // now we know the entries of ILDL(1,droptol) w.r.t. L
         	     
      
      
      // we just started this block
      if (j==startblock) {
#ifdef PRINT_INFO1
	 printf("start block with column %ld\n",j);
#endif
	 // clear idxL and use idxLj instead
	 // the fastest way to achieve this is a simple swap
	 pi=idxL;    idxL=idxLj;       idxLj=pi;
	 i=cntL;     cntL=cntLj;       cntLj=i;
	 pi=idxposL; idxposL=idxposLj; idxposLj=pi;
	 // reset additional fill-in
	 cnt_fill=0;
	 // we are at the last step, finish this block anyway
	 if (j==n-1) {
#ifdef PRINT_INFO1
	    printf("final step, quit with final blocksize[%ld]=%ld\n",*nblocks,endblock-startblock+1);
#endif
	    blocksize[(*nblocks)++]=endblock-startblock+1;
	 }
      }
      else { // measure the additional fill when merging 
	 // old fill
	 k=endblock-startblock+1;
	 old_fill=(2*cntL+k)*k;
	 // additional fill by step j
	 add_fill=2*cntLj+1;
	 // compute the amount of fill when not merging
	 old_fill+=add_fill;
	 // new fill when merging the blocks
	 // increase block size
	 k=k+1;
	 // augment block sub-diagonal and block super-diagonal patterns
	 // but exclude diagonal j entry from old fill
	 ii=0; jj=0;
	 l=0;
	 while (ii<cntL && jj<cntLj) {
	       t=idxL[ii]; s=idxLj[jj];
	       // make sure to exclude the diagonal entry
	       if (t==j)
		  ii++;
	       else {
		  if (t<s) {
		     buff[l++]=t;
		     ii++;
		  }
		  else if (s<t) {  
		     buff[l++]=s;
		     jj++;
		  }
		  else { // s=t
		     buff[l++]=t;
		     ii++;
		     jj++;
		  }
	       } // end if-else
	 } // end while
	 while (ii<cntL) {
	       t=idxL[ii];
	       // make sure to exclude the diagonal entry
	       if (t!=j)
		  buff[l++]=t;
	       ii++;
	 } // end while
	 while (jj<cntLj) {
	       s=idxLj[jj];
	       buff[l++]=s;
	       jj++;
	 } // end while
	 // make sure the j is removed from idxL
	 idxposL[j]=0;
	 // copy merged list back to idxL
	 cntL=0;
	 while (cntL<l) {
	       t=buff[cntL];
	       idxL[cntL]=t;
	       idxposL[t]=++cntL;
	 } // end while

	 new_fill=(2*cntL+k)*k;
	 // extend block
	 if (   new_fill<=ELBOW_BLOCK*(old_fill-cnt_fill)
	     || new_fill<=(old_fill-cnt_fill)+BLOCK_EXT*k
	     || (colj==j-1)) {
#ifdef PRINT_INFO1
	    printf("merge column %ld\n",j);
#endif
	    // shift end of the current block
	    endblock=endblock+1;
	    // increase the number of wasted memory
	    cnt_fill+=new_fill-old_fill;
	    // subsequent column was necessarily added
	    if (colj==j+1) {
	       // shift end of the current block again
	       endblock=endblock+1;
	       // increase the number of wasted memory
	       cnt_fill+=(2*cntL+(k+1))*(k+1)-(2*cntL+k)*k;
	       // advance beyond the subsequent column
	       j++;
	    }
	    // we are at the last step, finish this block anyway
	    if (j==n-1) {
#ifdef PRINT_INFO1
	       printf("final step, quit with final blocksize[%ld]=%ld\n",*nblocks,endblock-startblock+1);
#endif
	       blocksize[(*nblocks)++]=endblock-startblock+1;
	    }
	 }
	 else {// store current block and start new block
#ifdef PRINT_INFO1
	    printf("end block with column %ld, do not merge column %ld, use blocksize[%ld]=%ld\n",endblock,j,*nblocks,endblock-startblock+1);
#endif
	    blocksize[(*nblocks)++]=endblock-startblock+1;
	    // clear idxL, idxU and use idxLj, idxUj instead
	    // the fastest way to achieve this is a simple swap
	    pi=idxL;    idxL=idxLj;       idxLj=pi;
	    i=cntL;     cntL=cntLj;       cntLj=i;
	    pi=idxposL; idxposL=idxposLj; idxposLj=pi;
#ifdef PRINT_INFO1
	    printf("start block with column %ld\n",j);
#endif
	    startblock=j; endblock=j;
	    // reset fill
	    cnt_fill=0;
	    if (j==n-1) {
#ifdef PRINT_INFO1
	       printf("final step, quit with final blocksize[%ld]=%ld\n",*nblocks,endblock-startblock+1);
#endif
	       blocksize[(*nblocks)++]=endblock-startblock+1;
	    }
	 } // end if-else
      } //end if-else j=startblock
      
      // reset idxposLj
      for (m=0; m<cntLj; m++) {
	  t=idxLj[m];
	  idxposLj[t]=0;
      } // end for m
      cntLj=0;
  } // end for j
  // **************************************************************************
  // *****                     END main loop                              *****
  // **************************************************************************


  // post processing
  j=0;
  for (i=0; i<*nblocks; i++)
      j+=blocksize[i];
#ifdef PRINT_INFO
  // printf("subdiagonal part\n");
  // for (i=0; i<n-1; i++)
  //    printf("%12.4le",Asdiag[i]);
  // printf("\n");
  // fflush(stdout);
  if (SL!=NULL) {
     printf("scaling\n");
     for (i=0; i<n; i++)
         printf("%12.4le",SL[i]);
     printf("\n");
     fflush(stdout);
  }
  printf("permutation\n");
  for (i=0; i<n; i++)
      printf("%6ld",p[i]);
  printf("\n");
  fflush(stdout);
#endif
#ifdef PRINT_INFO
  printf("block partitioning, number of blocks: %ld, check sum %ld\n",*nblocks,j);
  for (i=0; i<*nblocks; i++)
      printf("%4ld",blocksize[i]);
  printf("\n");
  fflush(stdout);
#endif


  // release memory
  free(idxL);
  free(idxposL);
  free(idxLj);
  free(idxposLj);

  free(buff);
  
  free(Afirst);
  free(ATfirst);
  free(Adiag);
  free(Asdiag);

  for (i=0; i<n; i++) {
      free(AT.rowind[i]);
      free(AT.val[i]);
  }
  free(AT.ncol);
  free(AT.rowind);
  free(AT.val);
  
  return (n-j);
} // end bildl1t_blocks

 

 
