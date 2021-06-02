/* $Id: ilu1t_blocks.c 5106 2019-06-19 16:54:15Z bolle $ 

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

integer ILU1T_BLOCKS(SPARSEMATRIX *A,
		     REALS *SL, REALS *SR, integer *p, integer *invq,
		     integer *blocksize, integer *nblocks,
		     REALS droptol)
{  
  /*
    given a reordering and scaling determine approximately a block structure
    by locally simulating ILU(1,droptol), i.e., we only allow fill from the
    neighbouring original nodes and drop small entries in modulus


    Input
    -----
    A           assumed to be a nonsingular sparse matrix 
    p,invq      permutation p and inverse permutation invq associated with the a priori
                permutation matrices P and Q
		we further assume that the row indices of A(q,p) in each column k=p(i)
		are already sorted in increasing order. 
    SL,SR       real diagonal scaling matrices

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
          cntLj_old,    // old length of idxLj when needed

          *idxU,       // temporary index list of column indices
          *idxposU,    // associated inverse mapping storing the position
          cntU,        // length of idxU
          cntU_old,    // old length of idxU when needed

          *idxUj,       // temporary index list of column indices row j
          *idxposUj,    // associated inverse mapping storing the position
          cntUj,        // length of idxUj
          cntUj_old;    // old length of idxUj when needed        

  
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2;     // temporary numerical pointers

  REALS   ajj,att,atj,ajt,ajtjt,aitit, // temporary matrix entries by modulus
          *pr,*pr2,      // temporary real numerical pointers
          rval,         // temporary scalar real numerical value
          *Adiag;       // array of diagonal entries in modulus
  DSparseMatrix AT;
  

  
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposL=(integer *)calloc(n,sizeof(integer));
  idxU   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposU=(integer *)calloc(n,sizeof(integer));
  // create additional temporary index lists at step jfor rows and column
  idxLj   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposLj=(integer *)calloc(n,sizeof(integer));
  idxUj   =(integer *)malloc(n*sizeof(integer));
  // init with zeros
  idxposUj=(integer *)calloc(n,sizeof(integer));
  
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

  // prepare memory for A^T
  AT.nr=AT.nc=n;
  AT.ncol  =(integer *) calloc(n,sizeof(integer));
  AT.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
  AT.val   =(REALS **)  malloc((size_t)n*sizeof(REALS *));

  // count nnz for A^T
  pi=AT.ncol;
  for (i=0; i<n; i++) {
      // scan column A(q,p[i])
      // shortcuts
      ii=p[i]; 
      pi2=A->rowind[ii];
      l=A->ncol[ii];
      for (j=0; j<l; j++) {
	  // index k of A(q[k],p[i])
	  kk=pi2[j];
	  // q[k]=kk
	  k=invq[kk];
	  // increase number of elements in row k
	  (pi[k])++;
      } // end for j
  } // end for i

  // set up remaining memory for AT=A(q,p)^T
  for (i=0; i<n; i++) {
      // nz in row i of A(q,p)
      l=pi[i];
      AT.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
      AT.val[i]   =(REALS *)  malloc((size_t)l*sizeof(REALS));
      // reset AT.ncol
      pi[i]=0;
  } // end for i

  // copy |SL*A(q,p)*SR| to A^T
  if (SL==NULL && SR==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(q,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(q[k],p[i])
	     kk=pi2[j];
	     // q[k]=kk
	     k=invq[kk];
	     // |A(q[k],p[i])|
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
  else if (SR==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(q,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(q[k],p[i])
	     kk=pi2[j];
	     // q[k]=kk
	     k=invq[kk];
	     // SL(q[k],q[k])|A(q[k],p[i])|
	     rval=SL[kk]*FABS(pA[j]);
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     AT.val[k][m]=rval;
	     pi[k]=m+1;
	 } // end for j
     } // end for i
  }
  else if (SL==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(q,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(q[k],p[i])
	     kk=pi2[j];
	     // q[k]=kk
	     k=invq[kk];
	     // |A(q[k],p[i])|SR(p[i],p[i])
	     rval=FABS(pA[j])*SR[ii];
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
  else { // SL!=0 and SR!=0
     for (i=0; i<n; i++) {
         // scan column A(q,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(q[k],p[i])
	     kk=pi2[j];
	     // q[k]=kk
	     k=invq[kk];
	     // SL(q[k],q[k]) |A(q[k],p[i])| SR(p[i],p[i])
	     rval=SL[kk]*FABS(pA[j])*SR[ii];
	     // current position in AT(:,k)
	     m=pi[k];
	     // store index and value
	     AT.rowind[k][m]=i;
	     AT.val[k][m]=rval;
	     // increase nz in column k of AT
	     pi[k]=m+1;
	 } // end for j
     } // end for i
  } // end if-elseif-else SL=0 and SR=0

  // scan A(q,p), A(q,p)^T in order to find the position of the diagonal entry
  // though it exists (otherwise use position that follows the digonal entry)
  for (i=0; i<n; i++) {
      // scan column A(q,p[i])
      // shortcuts
      ii=p[i]; 
      pi2=A->rowind[ii];
      l=A->ncol[ii];
      for (j=0; j<l; j++) {
	  // index k of A(q[k],p[i])
	  kk=pi2[j];
	  // q[k]=kk
	  k=invq[kk];
	  // diagonal and sub-diagonal begins
	  if (k>=i) {
	     break;
	  } // end if k>=i
      } // end for j
      Afirst[i]=j;
      // scan row A(q[i],p)
      // shortcuts
      pi2=AT.rowind[i];
      pr=AT.val[i];
      l=AT.ncol[i];
      for (j=0; j<l; j++) {
	  // index k of A(q[i],p[k])
	  k=pi2[j];
	  // diagonal and super-diagonal begins
	  if (k>=i) {
	     // |A(q,p)| found
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
  cntL=cntU=0;   // counter for accumulated index lists
  cntLj=cntUj=0; // counter for temporary current index lists
  startblock=endblock=0; // initial beginning of the current block
  cnt_fill=0; // counter for additional fill-in when merging blocks
  *nblocks=0; // number of diagonal blocks of A
  for (j=0; j<n; j++) {

      // compute fill-in in column j of L, row j of U
	    
#ifdef PRINT_INFO
      printf("startblock=%3ld, endblock=%3ld\n",startblock,endblock);
      for (m=0; m<cntL; m++)
	  printf("%4ld",idxL[m]);
      printf("\n");
      for (m=0; m<cntU; m++)
	  printf("%4ld",idxU[m]);
      printf("\n");
      printf("column/row j=%3ld\n",j);
      fflush(stdout);
#endif

      // threshold*|A(q[j],p[j])|
      ajj=droptol*Adiag[j];
      
      // -----------------------------------------------
      // -----------------------------------------------
      // ------------ extract COLUMN j of A(q,p) -------
      // check sub-diagonal part of COLUMN j of A(q,p) for its numerical fill
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
	  // tt=q[t]
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
	     }
	     if (SR!=NULL) { // SL!=NULL
	        rval=SR[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		val*=rval;
#else
		val.r*=rval;
		val.i*=rval;
#endif
	     }

	     // check whether the entry is greater than a threshold relative to the diagonal entry
	     if (FABS(val)>=ajj) {
	        // a fill entry
	        idxLj[cntLj]=t;       // bookmark index and ...
		idxposLj[t]=++cntLj;  // its physical position, shifted by 1
	     } // end if 
	  } // end if t>j
      } // end for m
#ifdef PRINT_INFO1
      printf("L(:,%3ld): fill ILU(0,%8.1le)\n",j,droptol);
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      
      // now we know the entries of ILU(0,droptol) w.r.t. L
      // next proceed with ILU(1,droptol)
      // check COLUMNS t<j of A(q,p) for their numerical fill
      // the entries 0<=t<j are located before r=Afirst[j]
#ifdef PRINT_INFO1
      printf("Afirst[%ld]=%ld\n",j,Afirst[j]);
#endif
      for (m=0; m<r && cntLj<(n-j-1); m++) {
	  // row index tt of A(tt,p[j])
	  tt=pi[m];
	  // tt=q[t]
	  t=invq[tt];
#ifdef PRINT_INFO1
	  printf("t=%ld\n",t);
#endif

	  // store counter before entering column t<j
	  cntLj_old=cntLj;
	  // super-diagonal entry, include column t of A(q,p) for ILU(1,droptol)
	  if (t<j) {
#ifdef PRINT_INFO1
	  printf("... also check column t=%3ld\n",t);
	  fflush(stdout);
#endif
	     // extract associated numerical value A(q[t],p[j])
	     val=pA[m];
	     if (SL!=NULL) {
	        rval=SL[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		val*=rval;
#else
		val.r=rval;
		val.i=rval;
#endif
	     } // end if SL!=0
	     if (SR!=NULL) {
	        rval=SR[jj];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
	        val*=rval;
#else
		val.r*=rval;
		val.i*=rval;
#endif
	     } // end if SR!=0
	     // |A(q[t],p[j])|
	     atj=FABS(val);
	     
	     // threshold |A(q[t],p[t])A(q[j],p[j])|
	     att=Adiag[t]*ajj;

	     // search for beginning of A(q[j:end],p[t])
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
		 // ss=q[s]
		 s=invq[ss];
		 // trailing part
		 if (s>j)
		    break;
	     } // end for l
	     Afirst[t]=l;
#ifdef PRINT_INFO1
	     printf("Afirst[%ld]=%ld\n",t,Afirst[t]);
#endif
	     
	     // check COLUMN t<j of A(q,p) for fill-in
	     // right now, all t>j are located at position Afirst[t] and higher
	     for (; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         ss=pi2[l];
		 // ss=q[s]
		 s=invq[ss];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposLj[s]) {	     
		    // extract associated numerical value A(q[s],p[t])
		    val=pA2[l];
		    if (SL!=NULL) {
		       rval=SL[ss];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		    } // end if SL!=0
		    if (SR!=NULL) {
		       rval=SR[tt];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		
		       val*=rval;
#else
		       val.r*=rval;
		       val.i*=rval;
#endif
		    } // end if SR!=0

		    // threshold relative to the diagonal entry
		    // |A(q[s],p[t])A(q[t],p[j])| >= droptol |A(q[t],p[t])A(q[j],p[j])|
		    if (FABS(val)*atj>=att) {
		       // a fill entry
		       idxLj[cntLj]=s;       // bookmark index and ...
		       idxposLj[s]=++cntLj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l
	     // fill of column t incorporated
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
      printf("L(:,%3ld): fill ILU(1,%8.1le)\n",j,droptol);
      for (m=0; m<cntLj; m++)
	  printf("%8ld",idxLj[m]);
      printf("\n");
      fflush(stdout);
#endif
      // -----------------------------------------------
      // -----------------------------------------------
      // ---------- END extract COLUMN j of A ----------

      
      // now we know the entries of ILU(1,droptol) w.r.t. L
      // next to similar computations for U

      
      // -----------------------------------------------
      // -----------------------------------------------
      // ------------ extract ROW j of A ------------

      // check super-diagonal part of ROW j of A(q,p) for its numerical fill
      // ATfirst[j] starts with indices t>=j
      // shortcuts
      pi=AT.rowind[j];
      pr=AT.val[j];
      l=AT.ncol[j];
      r=ATfirst[j];
      for (m=r; m<l; m++) {
	  // column index t of A(q[j],p[t])^T
	  t=pi[m];

	  // super-diagonal fill
	  if (t>j) {
	     // extract associated numerical value
	     rval=pr[m];

	     // check whether the entry is greater than a threshold relative to the diagonal entry
	     if (rval>=ajj) {
	        // a fill entry
	        idxUj[cntUj]=t;       // bookmark index and ...
		idxposUj[t]=++cntUj;  // its physical position, shifted by 1
	     } // end if 
	  } // end if t>j
      } // end for m
#ifdef PRINT_INFO1
      printf("U(%3ld,:): fill ILU(0,%8.1le)\n",j,droptol);
      for (m=0; m<cntUj; m++)
	  printf("%8ld",idxUj[m]);
      printf("\n");
      fflush(stdout);
#endif
      
      // now we know the entries of ILU(0,droptol) w.r.t. U
      // next proceed with ILU(1,droptol)
      // check ROWS t<j of A(q,p) for their numerical fill
      // the entries 0<=t<j are located before t=ATfirst[j]
#ifdef PRINT_INFO1
      printf("ATfirst[%ld]=%ld\n",j,ATfirst[j]);
#endif
      for (m=0; m<r && cntUj<(n-j-1); m++) {
	  // column index t of A(q[j],p[t])
	  t=pi[m];

	  // store counter before entering row t<j
	  cntUj_old=cntUj;
	  // super-diagonal entry, include row t of A(q,p) for ILU(1,droptol)
	  if (t<j) {
#ifdef PRINT_INFO1
	     printf("... also check row t=%3ld\n",t);
	     fflush(stdout);
#endif
	     // extract associated numerical value |A(q[j],p[t])|
	     ajt=pr[m];
	     

	     // threshold |A(q[t],p[t])A(q[j],p[j])|
	     att=Adiag[t]*ajj;

	     // search for beginning of A(q[t],p[j:end])
	     // right now ATfirst[t] is only guaranteed to start with some
	     // index s>=t
	     // shortcuts
	     pi2=AT.rowind[t];
	     pr2=AT.val[t];
	     ii=AT.ncol[t];
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A(ss,p[t])
	         s=pi2[l];
		 // trailing part
		 if (s>j) {	     
		    break;
		 }
	     } // end for l
	     ATfirst[t]=l;
#ifdef PRINT_INFO1
	     printf("ATfirst[%ld]=%ld\n",t,ATfirst[t]);
#endif
	     
	     // check ROW t<j of A(q,p) for fill-in
	     // right now, all t>j are located at position ATfirst[t] and higher
	     for (l=ATfirst[t]; l<ii; l++) {
	         // row index ss of A(q[t],p[s])
	         s=pi2[l];
	  
		 // sub-diagonal entry that has not been considered as fill yet
		 if (!idxposUj[s]) {	     
		    // extract associated numerical value A(q[t],p[s])
		    rval=pr2[l];

		    // threshold relative to the diagonal entry
		    // |A(q[j],p[t])A(q[t],p[s])| >= droptol |A(q[t],p[t])A(q[j],p[j])|
		    if (rval*ajt>=att) {
		       // a fill entry
		       idxUj[cntUj]=s;       // bookmark index and ...
		       idxposUj[s]=++cntUj;  // its physical position, shifted by 1
		    } // end if 
		 } // end if s>j
	     } // end for l
	     // fill of row t incorporated
	     // now idxUj[0,...,cntUj_old-1] and idxUj[cntUj_old,...,cntUj-1]
	     // are distinct index arrays each with indices in increasing order
	     // we temporarily merge them in buff if necessary
	     if (0<cntUj_old && cntUj_old<cntUj) {
	        ii=0; jj=cntUj_old;
		l=0;
		while (ii<cntUj_old && jj<cntUj) {
	              t=idxUj[ii]; tt=idxUj[jj];
		      if (t<tt) {
			 buff[l++]=t;
			 ii++;
		      }
		      else { // tt<t, the case tt=t is avoided by flagging idxposLj
			 buff[l++]=tt;
			 jj++;
		      }
		} // end while
		while (ii<cntUj_old) {
		      t=idxUj[ii];
		      buff[l++]=t;
		      ii++;
		} // end while
		while (jj<cntUj) {
		      tt=idxUj[jj];
		      buff[l++]=tt;
		      jj++;
		} // end while
		// write the sorted list back to idxUj
		for (l=0; l<cntUj; ) {
		    t=buff[l];
		    idxUj[l]=t;
		    idxposUj[t]=++l;
		} // end for l
	     } // end if 0<cntUj_old<cntUj
	  } // end if t<j
      } // end for m
#ifdef PRINT_INFO
      printf("U(%3ld,:): fill ILU(1,%8.1le)\n",j,droptol);
      for (m=0; m<cntUj; m++)
	  printf("%8ld",idxUj[m]);
      printf("\n");
      fflush(stdout);
#endif
      // -----------------------------------------------
      // -----------------------------------------------
      // ---------- END extract ROW j of A ----------


      
      // we just started this block
      if (j==startblock) {
#ifdef PRINT_INFO1
	 printf("start block with column %ld\n",j);
#endif
	 // clear idxL, idxU and use idxLj, idxUj instead
	 // the fastest way to achieve this is a simple swap
	 pi=idxL;    idxL=idxLj;       idxLj=pi;
	 i=cntL;     cntL=cntLj;       cntLj=i;
	 pi=idxposL; idxposL=idxposLj; idxposLj=pi;
	 pi=idxU;    idxU=idxUj;       idxUj=pi;
	 i=cntU;     cntU=cntUj;       cntUj=i;
	 pi=idxposU; idxposU=idxposUj; idxposUj=pi;
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
	 old_fill=(cntL+cntU+k)*k;
	 // additional fill by step j
	 add_fill=cntLj+cntUj+1;
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

	 ii=0; jj=0;
	 l=0;
	 while (ii<cntU && jj<cntUj) {
	       t=idxU[ii]; s=idxUj[jj];
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
	 while (ii<cntU) {
	       t=idxU[ii];
	       // make sure to exclude the diagonal entry
	       if (t!=j)
		  buff[l++]=t;
	       ii++;
	 } // end while
	 while (jj<cntUj) {
	       s=idxUj[jj];
	       buff[l++]=s;
	       jj++;
	 } // end while
	 // make sure the j is removed from idxU
	 idxposU[j]=0;
	 // copy merged list back to idxU
	 cntU=0;
	 while (cntU<l) {
	       t=buff[cntU];
	       idxU[cntU]=t;
	       idxposU[t]=++cntU;
	 } // end while

	 new_fill=(cntL+cntU+k)*k;
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
	    pi=idxU;    idxU=idxUj;       idxUj=pi;
	    i=cntU;     cntU=cntUj;       cntUj=i;
	    pi=idxposU; idxposU=idxposUj; idxposUj=pi;
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
      // reset idxposUj
      for (m=0; m<cntUj; m++) {
	  t=idxUj[m];
	  idxposUj[t]=0;
      } // end for m
      cntUj=0;
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
  free(idxU);
  free(idxposU);
  free(idxLj);
  free(idxposLj);
  free(idxUj);
  free(idxposUj);

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

 

 
