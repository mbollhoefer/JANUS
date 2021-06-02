/* $Id: ildl1t_blocks.c 4145 2018-04-19 11:55:02Z bolle $ */
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


#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define ELBOW_BLOCK     1.3333
#define BLOCK_EXT       4


// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

integer ILDL1T_BLOCKS(SPARSEMATRIX *A,
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

          *idxLj,       // temporary index list of row indices column j
          *idxposLj,    // associated inverse mapping storing the position
          cntLj,        // length of idxLj
          cntLj_old;    // old length of idxLj when needed

  
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2;     // temporary numerical pointers

  REALS   ajj,att,atj,ajt,ajtjt,aitit, // temporary matrix entries by modulus
          *pr,*pr2, *pAT, // temporary real numerical pointers
          rval,         // temporary scalar real numerical value
          *Adiag;       // array of diagonal entries in modulus
  DSparseMatrix AT;
  

  
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposL=(integer *)calloc(n,sizeof(integer));
  // create additional temporary index lists at step jfor rows and column
  idxLj   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposLj=(integer *)calloc(n,sizeof(integer));
  
  // buffer for index lists
  nbuff =n;
  buff  =(integer *)malloc((size_t)nbuff *sizeof(integer));

  // beginning of each column/row at step j
  Afirst =(integer *)malloc((size_t)n*sizeof(integer));
  ATfirst=(integer *)malloc((size_t)n*sizeof(integer));
  // diagonal entries in modulus
  Adiag  =(REALS *)malloc(n*sizeof(REALS));
  for (i=0; i<n; i++)
      Adiag[i]=0.0;

  // prepare memory for A^T, recall that for symmetry reasons, only half of A is stored
  AT.nr=AT.nc=n;
  AT.ncol  =(integer *) calloc(n,sizeof(integer));
  AT.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
  AT.val   =(REALS **)  malloc((size_t)n*sizeof(REALS *));

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

  // set up remaining memory for AT=A(p,p)^T
  for (i=0; i<n; i++) {
      // nz in row i of A(p,p)
      l=pi[i];
      AT.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
      AT.val[i]   =(REALS *)  malloc((size_t)l*sizeof(REALS));
      // reset AT.ncol
      pi[i]=0;
  } // end for i

  // copy |SL*A(p,p)^T*SL| to A^T
  if (SL==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // |A(p[k],p[i])|
	     rval=FABS(pA[j]);
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     AT.val[k][m]=rval;
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
	 for (j=0; j<l; j++) {
	     // index k of A(p[k],p[i])
	     kk=pi2[j];
	     // p[k]=kk
	     k=invq[kk];
	     // SL(p[k],p[k]) |A(p[k],p[i])| SL(p[i],p[i])
	     rval=SL[kk]*FABS(pA[j])*SL[ii];
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     AT.val[k][m]=rval;
	     // increase nz in column k of AT
	     pi[k]=m+1;
	 } // end for j
     } // end for i
  } // end if-else SL=0

  // scan A(p,p), A(p,p)^T in order to find the position of the diagonal entry
  // though it exists (otherwise use position that follows the digonal entry)
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
	  // diagonal and sub-diagonal begins
	  if (k>=i) {
	     break;
	  } // end if k>=i
      } // end for j
      Afirst[i]=j;
      // scan row A(p[i],p)
      // shortcuts
      pi2=AT.rowind[i];
      pr=AT.val[i];
      l=AT.ncol[i];
      for (j=0; j<l; j++) {
	  // index k of A(p[i],p[k])
	  k=pi2[j];
	  // diagonal and super-diagonal begins
	  if (k>=i) {
	     // |A(p,p)| found
	     if (k==i)
	        Adiag[i]=pr[j];
	     // store the beginning of the diagonal/super-diagonal part
	     break;
	  } // end if k>=i
      } // end for j
      ATfirst[i]=j;
  } // end for i
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
  cntLj=0; // counter for temporary current index lists
  startblock=endblock=0; // initial beginning of the current block
  cnt_fill=0; // counter for additional fill-in when merging blocks
  *nblocks=0; // number of diagonal blocks of A
  for (j=0; j<n; j++) {

      // compute fill-in in column j of L
	    
#ifdef PRINT_INFO
      printf("startblock=%3ld, endblock=%3ld\n",startblock,endblock);
      for (m=0; m<cntL; m++)
	  printf("%4ld",idxL[m]);
      printf("\n");
      printf("column j=%3ld\n",j);
      fflush(stdout);
#endif

      // threshold*|A(p[j],p[j])|
      ajj=droptol*Adiag[j];
      
      // -----------------------------------------------
      // -----------------------------------------------
      // ------------ extract COLUMN j of A(q,p) -------
      // check sub-diagonal part of COLUMN j of A(p,p) for its numerical fill
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
	  if (t>j) {
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
	     if (FABS(val)>=ajj) {
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
      // ------------ extract COLUMN j of AT=|SL*A(p,p)*SL|^T ------------
      // check sub-diagonal part of COLUMN j of |SL*A(p,p)*SL|^T for its numerical fill
      // ATfirst[j] starts with indices t>=j
      // shortcuts
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
	  if (t>j) {
	     // check whether the entry is greater than a threshold relative to the diagonal entry
	     if (pAT[m]>=ajj) {
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
#ifdef PRINT_INFO1
      printf("L(:,%3ld): fill ILDL(0,%8.1le) including AT=|SL*A(p,p)*SL|^T\n",
	     j,droptol);
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      
      // now we know the entries of ILDL(0,droptol) w.r.t. L
      // next proceed with ILDL(1,droptol)
      // check COLUMNS t<j of A(p,p) for their numerical fill
      // the entries 0<=t<j are located before r=Afirst[j]
#ifdef PRINT_INFO1
      printf("Afirst[%ld]=%ld\n",j,Afirst[j]);
#endif
      jj=p[j]; 
      pi=A->rowind[jj];
      r=Afirst[j];
      for (m=0; m<r && cntLj<(n-j-1); m++) {
	  // row index tt of A(tt,p[j])
	  tt=pi[m];
	  // tt=p[t]
	  t=invq[tt];
#ifdef PRINT_INFO1
	  printf("t=%ld\n",t);
#endif

	  // super-diagonal entry, include column t of A(p,p) for ILDL(1,droptol)
	  if (t<j) {
#ifdef PRINT_INFO1
	     printf("... also check column t=%3ld\n",t);
	     fflush(stdout);
#endif
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
	     
	     // threshold |A(p[t],p[t])A(p[j],p[j])|
	     att=Adiag[t]*ajj;

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
		 if (s>j)
		    break;
	     } // end for l
	     Afirst[t]=l;
#ifdef PRINT_INFO1
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

		    // threshold relative to the diagonal entry
		    // |A(p[s],p[t])A(p[t],p[j])| >= droptol |A(p[t],p[t])A(p[j],p[j])|
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
	     pr2=AT.val[t];
	     ii=AT.ncol[t];
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A^T(p[s],p[t])
	         s=pi2[l];
		 // trailing part
		 if (s>j)
		    break;
	     } // end for l
	     ATfirst[t]=l;
#ifdef PRINT_INFO1
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
		    // threshold relative to the diagonal entry
		    // |A^T(p[s],p[t])A(p[t],p[j])| >= droptol |A(p[t],p[t])A(p[j],p[j])|
		    if (pr2[l]*atj>=att) {
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
	  } // end if t<j
      } // end for m

      // check COLUMNS t<j of A^T(p,p) for their numerical fill
      // the entries 0<=t<j are located before r=ATfirst[j]
#ifdef PRINT_INFO1
      printf("ATfirst[%ld]=%ld\n",j,ATfirst[j]);
#endif
      pi=AT.rowind[j];
      pAT=AT.val[j];
      r=ATfirst[j];
      for (m=0; m<r && cntLj<(n-j-1); m++) {
	  // row index t of A(p[t],p[j])
	  t=pi[m];
#ifdef PRINT_INFO1
	  printf("t=%ld\n",t);
#endif

	  // super-diagonal entry, include column t of A^T(p,p) for ILDL(1,droptol)
	  if (t<j) {
#ifdef PRINT_INFO1
	     printf("... also check column t=%3ld\n",t);
	     fflush(stdout);
#endif
	     // extract associated numerical value A^T(p[t],p[j])
	     // |A^T(p[t],p[j])|
	     atj=pAT[m];
	     
	     // threshold |A(p[t],p[t])A(p[j],p[j])|
	     att=Adiag[t]*ajj;

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
		 if (s>j)
		    break;
	     } // end for l
	     Afirst[t]=l;
#ifdef PRINT_INFO1
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

		    // threshold relative to the diagonal entry
		    // |A(p[s],p[t])A(p[t],p[j])| >= droptol |A(p[t],p[t])A(p[j],p[j])|
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
	     pr2=AT.val[t];
	     ii=AT.ncol[t];
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A^T(p[s],p[t])
	         s=pi2[l];
		 // trailing part
		 if (s>j)
		    break;
	     } // end for l
	     ATfirst[t]=l;
#ifdef PRINT_INFO1
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
		    // threshold relative to the diagonal entry
		    // |A^T(p[s],p[t])A(p[t],p[j])| >= droptol |A(p[t],p[t])A(p[j],p[j])|
		    if (pr2[l]*atj>=att) {
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
	  } // end if t<j
      } // end for m

      
#ifdef PRINT_INFO
      printf("L(:,%3ld): fill ILDL(1,%8.1le)\n",j,droptol);
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
	     || new_fill<=(old_fill-cnt_fill)+BLOCK_EXT*k) {
#ifdef PRINT_INFO1
	    printf("merge column %ld\n",j);
#endif
	    // shift end of the current block
	    endblock=endblock+1;
	    // increase the number of wasted memory
	    cnt_fill+=new_fill-old_fill;
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
  
  // remove tiny blocks, this increases the number of blocks
  i=0;
  for (j=0; j<*nblocks; j++) {
      l=blocksize[j];
      if (l<MIN_BLOCK_SIZE_JANUS)
	 i+=(l-1);
  }
  // new final position in blocksize
  i+=*nblocks-1;
  for (j=*nblocks-1; j>=0; j--) {
      l=blocksize[j];
      if (l<MIN_BLOCK_SIZE_JANUS) {
	 *nblocks+=l-1;
	 // use blocksize[j] times blocks of size 1
	 for (k=l; k>0; k--) 
	     blocksize[i--]=1;
      }
      else // simply shift
	 blocksize[i--]=l;
  } // end for j
  
  j=0;
  for (i=0; i<*nblocks; i++)
      j+=blocksize[i];
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

  for (i=0; i<n; i++) {
      free(AT.rowind[i]);
      free(AT.val[i]);
  }
  free(AT.ncol);
  free(AT.rowind);
  free(AT.val);
  
  return (n-j);
} // end ilu1t_blocks

 

 
