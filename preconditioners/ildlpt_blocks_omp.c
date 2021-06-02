/* $Id: ildlpt_blocks_omp.c 6094 2020-03-28 16:17:50Z bolle $ 

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	June 19, 2019. JANUS Block ILU R1.1.  

    Notice:

	Copyright (c) 2019 by TU Braunschweig.  All Rights Reserved.

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
#include <omp.h>

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
#define MIN_COLUMN_SIZE 10


// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

integer ILDLPT_BLOCKS(SPARSEMATRIX *A,
		      REALS *SL, REALS *SR, integer *p, integer *invq,
		      integer *blocksize, integer *nblocks,
		      REALS droptol, integer level_of_fill)
{  
  /*
    given a reordering and scaling determine approximately a block structure
    by locally simulating ILU(level_of_fill,droptol), i.e., we only allow fill from the
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
    
          *idxL,       // temporary index list of row indices
          cntL,        // length of idxL
          cntL_old,    // old length of idxL when needed

          *idxLj,       // temporary index list of row indices column j
          cntLj,        // length of idxLj
          cntLj_old,    // old length of idxLj when needed

          **adjL,       // adjacency list for each column of L
          *nzL,         // number of nonzeros in the adjacency list for each column of L

          *idxL_old,    // array containing the sub block-diagonal indices of the previous block-column
          *idxL_current,// array containing the sub block-diagonal indices of the current block-column

          *ibuff;       // index array for check marks
    
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2;     // temporary numerical pointers

  REALS   ajj,att,atj,ajt,ajtjt,aitit, // temporary matrix entries by modulus
          *pr,*pr2,     // temporary real numerical pointers
          rval,         // temporary scalar real numerical value
          *Adiag;       // array of diagonal entries in modulus
  DSparseMatrix AT;

  integer *pathlen, *maxpathlen, *flag, *flaglist, *queue, nqueue, startqueue, endqueue, nflagged,nzj0,nzj;
  REALS  *factorn, *factord, weightn, weightd, fkn, fkd, tau, droptol_end;

  integer *pathlen_global, *flag_global, *flaglist_global, *queue_global,
          *idxLj_global, *idxUj_global;
  REALS  *factorn_global, *factord_global;

  integer max_threads=omp_get_max_threads(), mythreadnum;


  // statistical array, count the path length
  maxpathlen=(integer *)calloc((size_t)n,sizeof(integer));
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc((size_t)n*sizeof(integer));
  // create additional temporary index lists at step jfor rows and column
  idxLj_global=(integer *)malloc((size_t)n*max_threads*sizeof(integer));

  // adjacency list for each column of L
  adjL=(integer **)malloc((size_t)n*sizeof(integer *));
  // number of nonzeros in the adjacency list for each column of L
  nzL =(integer *) calloc((size_t)n,sizeof(integer));
  for (j=0; j<n; j++) {
      adjL[j]=(integer *)malloc((size_t)MIN_COLUMN_SIZE*sizeof(integer));
  } // end for i
  
  // buffer for index lists
  nbuff =n;
  buff  =(integer *)malloc((size_t)nbuff *sizeof(integer));

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

  // copy |SL*A(p,p)*SL| to A^T
  if (SL==NULL) {
     for (i=0; i<n; i++) {
         // scan column A(p,p[i])
         // shortcuts
         ii=p[i]; 
	 pi2=A->rowind[ii];
	 pA=A->val[ii];
	 l=A->ncol[ii];
	 for (j=0; j<l; j++) {
	     // index k of A(q[k],p[i])
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
  } // end if-elseif-else SL=0

  // scan A(p,p)^T in order to find the position of the diagonal entry
  // though it exists 
  for (i=0; i<n; i++) {
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
  } // end for i

  pathlen_global =(integer *)calloc((size_t)n*max_threads,sizeof(integer));
  flag_global    =(integer *)calloc((size_t)n*max_threads,sizeof(integer));
  flaglist_global=(integer *)malloc((size_t)n*max_threads*sizeof(integer));
  queue_global   =(integer *)malloc((size_t)n*max_threads*sizeof(integer));
  factorn_global =(REALS *)  malloc((size_t)n*max_threads*sizeof(REALS));
  factord_global =(REALS *)  malloc((size_t)n*max_threads*sizeof(REALS));
  for (i=0; i<n*max_threads; i++)
      factorn_global[i]=factord_global[i]=1.0;

#ifdef PRINT_INFO1
  printf("(SL*A*SL)(p,p)\n");
  for (j=0; j<n; j++) {
       jj=p[j]; 
       pi=A->rowind[jj];
       pA=A->val[jj];
       m=A->ncol[jj];
       printf("column %5ld\n",j+1);
       for (r=0; r<m; r++) {
	   // row index tt of A(tt,p[k])
	   tt=pi[r];
	   // tt=p[l]
	   l=invq[tt];
	   printf("%8ld",l+1);
       } // end for r
       printf("\n"); fflush(stdout);
       for (r=0; r<m; r++) {
	   printf("%8.1le",pA[r]);
       } // end for r
       printf("\n"); fflush(stdout);
  }
  printf("\n");
#endif

  // droptol_end=droptol/log_2(n)
  droptol_end=droptol*log(2.0)/log((double)n);
  // **************************************************************************
  // *****                       main loop                                *****
  // **************************************************************************
#pragma omp parallel for\
  shared (n,idxLj_global,pathlen_global,flag_global,maxpathlen,\
	  flaglist_global,queue_global,factorn_global,factord_global,p,A,invq,\
	  SL,Adiag,droptol,droptol_end,level_of_fill,adjL,nzL,AT)\
  private (mythreadnum,idxLj,pathlen,flag,flaglist,queue,factorn,factord, \
	   nqueue,startqueue,endqueue,nflagged,cntLj,k,jj,pi,pA,m,r,tt,l,val,\
	   rval,weightn,weightd,pr,fkd,fkn,tau,nzj0,nzj)
  for (j=0; j<n; j++) {

      mythreadnum=omp_get_thread_num();

      idxLj=idxLj_global+n*mythreadnum;

      pathlen =pathlen_global +n*mythreadnum;
      flag    =flag_global    +n*mythreadnum;
      flaglist=flaglist_global+n*mythreadnum;
      queue   =queue_global   +n*mythreadnum;
      factorn =factorn_global +n*mythreadnum;
      factord =factord_global +n*mythreadnum;
    
      tau=(droptol*(n-1-j))/n+(droptol_end*j)/n;
      // initial amount of level-0 fill
      nzj0=1;
      // initial total amount of fill, weighted by the level
      nzj =1;
    
      // init queue with j
      queue[0]  =j;
      factord[0]=1.0;
      factorn[0]=1.0;
      nqueue=1;
      startqueue=0;
      endqueue=0;
      // flag j as visited
      flag[j]=1;
      nflagged=1;
      flaglist[0]=j;
      // number of adjacent nodes
      cntLj=0;      
      while (nqueue>0) {
	    // dequeue leading entry
	    k  =queue[startqueue];
	    fkd=factord[startqueue];
	    fkn=factorn[startqueue];
	    startqueue++; nqueue--; if (startqueue>=n) startqueue=0;
	    // search through all unvisited neighbours of k

	    // ------------ extract COLUMN k of tril(A(q,p)) -------
	    // recall that only half of the matrix is stored!
	    // shortcuts
	    jj=p[k]; 
	    pi=A->rowind[jj];
	    pA=A->val[jj];
	    m=A->ncol[jj];
	    for (r=0; r<m; r++) {
	        // row index tt of A(tt,p[k])
	        tt=pi[r];
		// tt=q[l]
		l=invq[tt];
		
		if (!flag[l]) {
		   val=pA[r];
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
		   }
		   // weight=fk * |A(p[l],p[k])/A(p[k],p[k])|
		   weightn=fkn*FABS(val);
		   weightd=fkd*Adiag[k];
		   if (weightn>=tau*(nzj0/(double)nzj)*weightd) {
		      // if (weightn>=tau*(nzj/(cntLj+1.0+pathlen[k]))*weightd) {
		      flag[l]=1;
		      flaglist[nflagged++]=l;
		      if (l<j && pathlen[k]<level_of_fill) {
			 // put l to the end of the queue
			 endqueue++; nqueue++; if (endqueue>=n) endqueue=0;
			 queue[endqueue]  =l;
			 factorn[endqueue]=weightn;
			 factord[endqueue]=weightd;
			 pathlen[l]=pathlen[k]+1;
			 maxpathlen[j]=MAX(maxpathlen[j],pathlen[l]);
		      } // end if
		      else if (l>j) {
			 idxLj[cntLj++]=l;
			 // amount of level-0 fill
			 if (pathlen[k]==0)
			    nzj0++;
			    // nzj++;
			 // total amount of fill, weighted by the level
			 nzj+=pathlen[k]+1;
		      } // end if-else if
		   } // end if weight>=tau
		} // end if !flag[l]
	    } // end for r
	    // ---------- END extract COLUMN k of tril(A(q,p)) -----

	    
	    // ----------- extract COLUMN k of tril(A(p,p))^T ------
	    // recall that only half of the matrix is stored!
	    // shortcuts
	    pi=AT.rowind[k];
	    pr=AT.val[k];
	    m=AT.ncol[k];
	    for (r=0; r<m; r++) {
	        // column index l of A(p[k],p[l])^T
	        l=pi[r];
		
		if (!flag[l]) {
		   // weight=fk * |A(p[k],p[l])^T/A(p[k],p[k])|
		   weightn=fkn*pr[r];
		   weightd=fkd*Adiag[k];
		   if (weightn>=tau*(nzj0/(double)nzj)*weightd) {
		      // if (weightn>=tau*(nzj/(cntLj+1.0+pathlen[k]))*weightd) {
		      flag[l]=1;
		      flaglist[nflagged++]=l;
		      if (l<j && pathlen[k]<level_of_fill) {
			 // put l to the end of the queue
			 endqueue++; nqueue++; if (endqueue>=n) endqueue=0;
			 queue[endqueue]  =l;
			 factorn[endqueue]=weightn;
			 factord[endqueue]=weightd;
			 pathlen[l]=pathlen[k]+1;
			 maxpathlen[j]=MAX(maxpathlen[j],pathlen[l]);
		      } // end if
		      else if (l>j) {
			 idxLj[cntLj++]=l;
			 // amount of level-0 fill
			 if (pathlen[k]==0)
			    nzj0++;
			    // nzj++;
			 // total amount of fill, weighted by the level
			 nzj+=pathlen[k]+1;
		      } // end if-else if
		   } // end if weight>=tau
		} // end if !flag[l]
	    } // end for r
	    // --------- END extract COLUMN k of tril(A(p,p))^T ----

      } // end while
      
      // sort entries in increasing order
      QQSORTI(idxLj,queue,&cntLj);
      // export idxLj to adjL[j]
      if (cntLj>MIN_COLUMN_SIZE)
	 adjL[j]=(integer *)realloc(adjL[j],(size_t)cntLj*sizeof(integer));
      memcpy(adjL[j],idxLj,(size_t)cntLj*sizeof(integer));
      nzL[j]=cntLj;
      // reset auxiliary arrays
      for (k=0; k<nflagged; k++) {
          l=flaglist[k];
	  flag[l]=0;
	  pathlen[l]=0;
      } // end  for j
      flag[j]=0;  
      pathlen[j]=0;

#ifdef PRINT_INFO1
      printf("adjL[%5ld]\n",j+1);
      for (k=0; k<nzL[j]; k++)
	  printf("%8ld",adjL[j][k]+1);
      printf("\n"); fflush(stdout);
#endif

  } // end for j


  
  free(idxLj_global);
  free(pathlen_global); 
  free(flag_global);
  free(flaglist_global);
  free(queue_global);
  free(factorn_global);
  free(factord_global);

  
  
  startblock=endblock=0; // initial beginning of the current block
  cnt_fill=0; // counter for additional fill-in when merging blocks
  *nblocks=0; // number of diagonal blocks of A
  cntL=0;   // counter for accumulated index lists
  for (j=0; j<n; j++) {
      
#ifdef PRINT_INFO1
      printf("startblock=%6ld, endblock=%6ld\n",startblock,endblock);
      for (m=0; m<cntL; m++)
	  printf("%6ld",idxL[m]);
      printf("\n");
      printf("current column/row j=%6ld\n",j);
      fflush(stdout);
#endif

      idxLj=adjL[j];
      cntLj=nzL[j];
      
      // we just started this block
      if (j==startblock) {
#ifdef PRINT_INFO1
	 printf("start block with column %ld\n",j);
#endif
	 // init idxL
	 cntL=cntLj;
	 memcpy(idxL,idxLj,(size_t)cntLj*sizeof(integer));
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
	 old_fill=(cntL+cntL+k)*k;
	 // additional fill by step j
	 add_fill=cntLj+cntLj+1;
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
	 // swap merged list with idxL
	 cntL=l;
	 pi=idxL; idxL=buff; buff=pi;


	 new_fill=(cntL+cntL+k)*k;
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

	    // init idxL
	    cntL=cntLj;
	    memcpy(idxL,idxLj,(size_t)cntLj*sizeof(integer));
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
      
  } // end for j
  // **************************************************************************
  // *****                     END main loop                              *****
  // **************************************************************************
#ifdef PRINT_INFO
#ifdef PRINT_INFO1
  j=0;
  for (i=0; i<*nblocks; i++)
      j+=blocksize[i];
  printf("initial block partitioning, number of blocks: %ld, check sum %ld\n",*nblocks,j);
  for (i=0; i<*nblocks; i++)
      printf("%4ld",blocksize[i]);
  printf("\n");
  fflush(stdout);
#endif
  r=0;
  rval=0.0;
  weightn=0.0;
  for (i=0; i<*nblocks; i++) {
      k=blocksize[i];

      r=MAX(r,k);
      rval+=k;
      weightn+=k*k;
  } /* end for i */

  // (*nblocks): number of blocks
  // r: maximum block size
  // rval: average block size
  rval=rval/(*nblocks);
  // weightd: standard deviation
  weightd=sqrt((weightn-*nblocks*rval*rval)/(*nblocks-1));
  printf("initial number of blocks: %ld, maximum block size: %ld, average %4.1lf(+-%8.1le)\n",*nblocks,r,rval,weightd);
#endif



  // now simulate progressive aggregation
  // auxiliary memory
  // array of indices for block-column j-1, sub block-diagonal part
  idxL_old    =(integer *)malloc((size_t)n*sizeof(integer));
  // array of indices for block-column j, sub block-diagonal part
  idxL_current=(integer *)malloc((size_t)n*sizeof(integer));
  // index buffer, must be initialized and re-set to 0
  ibuff       =(integer *)calloc((size_t)n,sizeof(integer));

  // counter for the current scalar row/column
  i=0;
  // counter for additional fill-in
  cnt_fill=0;
  // counters for the sub/super block-diagonal fill-in of block column j-1
  cntL_old=0;
  for (j=0; j<*nblocks; j++) {
      // current size of the diagonal block
      l=blocksize[j];
      // counters for the current sub/super block-diagonal fill-in, block-column j
      cntL=0;
      // merge fill inside block j
      for (k=i; k<i+l; k++) {
	  // --- sub block-diagonal part ---
	  // adjacency structure of column k
          idxLj=adjL[k];
	  // number of nonzeros in column k
          cntLj=nzL[k];
	  // merge sub block-diagonal entries of columns i,...,i+l-1 into common index array idxL_current
	  r=s=t=0;
	  while (r<cntL && s<cntLj)  {
	        // reduce indices to the sub block-diagonal part and merge
	        ii=idxL_current[r]; jj=idxLj[s];
		if (jj<i+l)
		   s++;
		else { // now ii,jj>=i+l, we have excluded the diagonal block
		       // by incremental construction, ii always fulfills ii>=i+l
		   if (ii<jj) {
		      buff[t++]=ii;
		      r++;
		   }
		   else if (jj<ii) {
		      buff[t++]=jj;
		      s++;
		   }
		   else { // jj=ii
		      buff[t++]=jj;
		      r++; s++;
		   }
		} // end if-else
	  } // end while
	  // copy remaining parts 
	  while (r<cntL)  {
	        // reduce I to its sub block-diagonal part
	        ii=idxL_current[r]; 
		// by incremental construction, ii always fulfills ii>=i+l
		buff[t++]=ii;
		r++;
	  } // end while
	  while (s<cntLj)  {
	        jj=idxLj[s];
		if (jj>=i+l) 
		   buff[t++]=jj;
		s++;
	  } // end while
	  // swap pointers buff, idxL_current
	  pi=buff; buff=idxL_current; idxL_current=pi;
	  cntL=t;
      } // end for k
      // now we have gathered the off-block-diagonal indices in block column/row j
      // into idxL_current with length cntL
      // now we can measure if merging block columns/rows j-1,j pays off
      
      // increase i by diagonal block size, i.e., i will advance to the next row/column behind block j
      i+=l;
      // check whether it makes sense to merge block j and j-1
      if (j>0) {
	 // first compute old fill, block j-1
	 k=blocksize[j-1];
	 ii=k*(k+2*cntL_old);
	 // l=blocksize[j]; // still valid
	 // add fill of block j
	 ii+=l*(l+2*cntL);
	 // --- block columns j-1 and j ---
	 // remove entries from block j-1 that will become part of the extended diagonal block j
	 // first step: check mark current nonzero entries
	 for (r=0; r<cntL; r++) {
	     s=idxL_current[r];
	     ibuff[s]=r+1;
	 } // end for r
	 // second step: compute additional entries from block j-1
	 // store them in buff, counter is t
	 t=0;
	 for (r=0; r<cntL_old; r++) {
	     s=idxL_old[r];
	     // scratch out entries from the diagonal block or those which are check-marked
	     if (s>=i && ibuff[s]==0)
	        buff[t++]=s;
	 } // end for r
	 // third step: clear check marks
	 for (r=0; r<cntL; r++) {
	     s=idxL_current[r];
	     ibuff[s]=0;
	 } // end for r
	 
	 // now we are able to compute the new fill
	 jj=(k+l)*(k+l+2*(cntL+t));
	 if ((jj<=ELBOW_BLOCK*(ii-cnt_fill) || jj<=(ii-cnt_fill)+BLOCK_EXT*MAX(k,l)) 
	     && k+l<=MAX_BLOCK_SIZE_JANUS) {
	    // incorporate additional fill-in
	    cnt_fill+=jj-ii;
	    blocksize[j]=k+l;
	    // delete block column j-1
	    blocksize[j-1]=0;
	    // merge block-column j-1 into block column j
	    // temporarily use idxL_old, since it is not needed anymore
	    r=s=kk=0;
	    while (r<cntL && s<t)  {
	          ii=idxL_current[r]; jj=buff[s];
		  if (ii<jj) {
		     idxL_old[kk++]=ii;
		     r++;
		  }
		  else { // jj>ii because of the check marks
		     idxL_old[kk++]=jj;
		     s++;
		  }
	    } // end while
	    while (r<cntL)  {
	          ii=idxL_current[r];
		  idxL_old[kk++]=ii;
		  r++;
	    } // end while
	    while (s<t)  {
	          jj=buff[s];
		  idxL_old[kk++]=jj;
		  s++;
	    } // end while
	    // swap pointers idxL_old, idxL_current
	    pi=idxL_old; idxL_old=idxL_current; idxL_current=pi;
	    cntL+=t;
	 }
	 else // reset fill counter
	    cnt_fill=0;
      } // end if j>0
      // shift blocks
      pi=idxL_current; idxL_current=idxL_old; idxL_old=pi;
      cntL_old=cntL;
  } // end for j


  // remove empty blocks
  j=0;
  i=0;
  k=0;
  for (j=0; j<*nblocks; j++) 
      if (blocksize[j]==0)
         k++;
      else
         blocksize[i++]=blocksize[j]; 
  *nblocks-=k;  

#ifdef PRINT_INFO
#ifdef PRINT_INFO1
  j=0;
  for (i=0; i<*nblocks; i++)
      j+=blocksize[i];
  printf("block partitioning after simulated progressive aggregation, number of blocks: %ld, check sum %ld\n",*nblocks,j);
  for (i=0; i<*nblocks; i++)
      printf("%4ld",blocksize[i]);
  printf("\n");
  fflush(stdout);
#endif
  r=0;
  rval=0.0;
  weightn=0.0;
  for (i=0; i<*nblocks; i++) {
      k=blocksize[i];

      r=MAX(r,k);
      rval+=k;
      weightn+=k*k;
  } /* end for i */

  // (*nblocks): number of blocks
  // r: maximum block size
  // rval: average block size
  rval=rval/(*nblocks);
  // weightd: standard deviation
  weightd=sqrt((weightn-*nblocks*rval*rval)/(*nblocks-1));
  printf("number of blocks after simulated progressive aggregation: %ld, maximum block size: %ld, average %4.1lf(+-%8.1le)\n",*nblocks,r,rval,weightd);
#endif

  
  
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
  
  //#ifdef PRINT_INFO
#ifdef PRINT_INFO1
  printf("block partitioning, number of blocks: %ld, check sum %ld\n",*nblocks,j);
  for (i=0; i<*nblocks; i++)
      printf("%4ld",blocksize[i]);
  printf("\n");
  fflush(stdout);
#endif
  r=0;
  rval=0.0;
  weightn=0.0;
  for (i=0; i<*nblocks; i++) {
      k=blocksize[i];

      r=MAX(r,k);
      rval+=k;
      weightn+=k*k;
  } /* end for i */

  // (*nblocks): number of blocks
  // r: maximum block size
  // rval: average block size
  rval=rval/(*nblocks);
  // weightd: standard deviation
  weightd=sqrt((weightn-*nblocks*rval*rval)/MAX(1,*nblocks-1));
  printf("number of blocks: %ld, maximum block size: %ld, average %4.1lf(+-%8.1le)\n",*nblocks,r,rval,weightd);

  r=0;
  rval=0.0;
  weightn=0.0;
  for (i=0; i<n; i++) {
      k=maxpathlen[i];

      r=MAX(r,k);
      rval+=k;
      weightn+=k*k;
  } /* end for i */

  // n: number of columns
  // r: maximum path length 
  // rval: average path length
  rval=rval/n;
  // weightd: standard deviation
  weightd=sqrt((weightn-*nblocks*rval*rval)/MAX(1,n-1));
  printf("size %ld, drop tolerance: %8.1le, levels of fill, maximum path length: %ld, average %4.1lf(+-%8.1le)\n",n,droptol,r,rval,weightd);
  //#endif


  // release memory
  free(idxL);

  free(buff);
  
  free(Adiag);

  for (i=0; i<n; i++) {
      free(AT.rowind[i]);
      free(AT.val[i]);
  }
  free(AT.ncol);
  free(AT.rowind);
  free(AT.val);

  for (i=0; i<n; i++) {
      if (adjL[i]!=NULL)
	 free(adjL[i]);
  } // end for i
  free(adjL);
  free(nzL);

  free(idxL_old);
  free(idxL_current);
  free(ibuff);
  free(maxpathlen);
  
  return (n-j);
} // end ilupt_blocks

 

 
