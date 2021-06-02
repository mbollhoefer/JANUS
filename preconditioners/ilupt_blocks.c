/* $Id: ilupt_blocks.c 5204 2019-07-17 16:45:00Z bolle $ 

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

integer ILUPT_BLOCKS(SPARSEMATRIX *A,
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

          *idxU,       // temporary index list of column indices
          cntU,        // length of idxU
          cntU_old,    // old length of idxU when needed

          *idxUj,       // temporary index list of column indices row j
          cntUj,        // length of idxUj
          cntUj_old,    // old length of idxUj when needed        
          **adjL,       // adjacency list for each column of L
          **adjU,       // adjacency list for each column of U^T    
          *nzL,         // number of nonzeros in the adjacency list for each column of L
          *nzU;         // number of nonzeros in the adjacency list for each column of U^T
    
  FLOAT   val,          // temporary scalar numerical value
          *pA,*pA2;     // temporary numerical pointers

  REALS   ajj,att,atj,ajt,ajtjt,aitit, // temporary matrix entries by modulus
          *pr,*pr2,     // temporary real numerical pointers
          rval,         // temporary scalar real numerical value
          *Adiag;       // array of diagonal entries in modulus
  DSparseMatrix AT;

  integer *pathlen, *flag, *flaglist, *queue, nqueue, startqueue, endqueue, nflagged,nzj0,nzj;
  double  *factorn, *factord, weightn, weightd, fkn, fkd, tau, droptol_end;
  

  
  // create temporary index lists for rows and column
  idxL   =(integer *)malloc(n*sizeof(integer));
  idxU   =(integer *)malloc(n*sizeof(integer));
  // create additional temporary index lists at step jfor rows and column
  idxLj   =(integer *)malloc(n*sizeof(integer));
  idxUj   =(integer *)malloc(n*sizeof(integer));

  // adjacency list for each column of L
  adjL=(integer **)malloc((size_t)n*sizeof(integer *));
  // adjacency list for each column of U^T    
  adjU=(integer **)malloc((size_t)n*sizeof(integer *));
  // number of nonzeros in the adjacency list for each column of L
  nzL =(integer *) calloc((size_t)n,sizeof(integer));
  // number of nonzeros in the adjacency list for each column of U^T
  nzU =(integer *) calloc((size_t)n,sizeof(integer));
  for (i=0; i<n; i++) {
      adjL[j]=(integer *)malloc((size_t)MIN_COLUMN_SIZE*sizeof(integer));
      adjU[j]=(integer *)malloc((size_t)MIN_COLUMN_SIZE*sizeof(integer));
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

  // scan A(q,p)^T in order to find the position of the diagonal entry
  // though it exists 
  for (i=0; i<n; i++) {
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
  } // end for i

  pathlen =(integer *)calloc((size_t)n,sizeof(integer));
  flag    =(integer *)calloc((size_t)n,sizeof(integer));
  flaglist=(integer *)malloc((size_t)n*sizeof(integer));
  queue   =(integer *)malloc((size_t)n*sizeof(integer));
  factorn =(double *) malloc((size_t)n*sizeof(double));
  factord =(double *) malloc((size_t)n*sizeof(double));
  for (i=0; i<n; i++)
      factorn[i]=factord[i]=1.0;

#ifdef PRINT_INFO1
  printf("(SL*A*SR)(q,p)\n");
  for (j=0; j<n; j++) {
       jj=p[j]; 
       pi=A->rowind[jj];
       pA=A->val[jj];
       m=A->ncol[jj];
       printf("column %5ld\n",j+1);
       for (r=0; r<m; r++) {
	   // row index tt of A(tt,p[k])
	   tt=pi[r];
	   // tt=q[l]
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
  for (j=0; j<n; j++) {
    
      tau=(droptol*(n-1-j))/n+(droptol_end*j)/n;
      // -----------------------------
      // --- lower triangular part ---
      // -----------------------------
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
	    // ------------ extract COLUMN k of A(q,p) -------
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
		   // weight=fk * |A(q[l],q[k])/A(q[k],p[k])|
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

      
    
    
      // -----------------------------
      // --- upper triangular part ---
      // -----------------------------
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
      cntUj=0;      
      while (nqueue>0) {
	    // dequeue leading entry
	    k  =queue[startqueue];
	    fkd=factord[startqueue];
	    fkn=factorn[startqueue];
	    startqueue++; nqueue--; if (startqueue>=n) startqueue=0;
	    // search through all unvisited neighbours of k
	    // ----------- extract COLUMN k of A(q,p)^T ------
	    // shortcuts
	    pi=AT.rowind[k];
	    pr=AT.val[k];
	    m=AT.ncol[k];
	    for (r=0; r<m; r++) {
	        // column index l of A(q[k],p[l])^T
	        l=pi[r];
		
		if (!flag[l]) {
		   // weight=fk * |A(q[k],p[l])^T/A(q[k],p[k])|
		   weightn=fkn*pr[r];
		   weightd=fkd*Adiag[k];
		   if (weightn>=tau*(nzj0/(double)nzj)*weightd) {
		      // if (weightn>=tau*(nzj/(cntUj+1.0+pathlen[k]))*weightd) {
		      flag[l]=1;
		      flaglist[nflagged++]=l;
		      if (l<j && pathlen[k]<level_of_fill) {
			 // put l to the end of the queue
			 endqueue++; nqueue++; if (endqueue>=n) endqueue=0;
			 queue[endqueue]  =l;
			 factorn[endqueue]=weightn;
			 factord[endqueue]=weightd;
			 pathlen[l]=pathlen[k]+1;
		      } // end if
		      else if (l>j) {
			 idxUj[cntUj++]=l;
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
      } // end while
      // sort entries in increasing order
      QQSORTI(idxUj,queue,&cntUj);
      // export idxLj to adjU[j]
      if (cntUj>MIN_COLUMN_SIZE)
	 adjU[j]=(integer *)realloc(adjU[j],(size_t)cntUj*sizeof(integer));
      memcpy(adjU[j],idxUj,(size_t)cntUj*sizeof(integer));
      nzU[j]=cntUj;
      // reset auxiliary arrays
      for (k=0; k<nflagged; k++) {
          l=flaglist[k];
	  flag[l]=0;
	  pathlen[l]=0;
      } // end  for j
      flag[j]=0;  
      pathlen[j]=0;

#ifdef PRINT_INFO1
      printf("adjU[%5ld]\n",j+1);
      for (k=0; k<nzU[j]; k++)
	  printf("%8ld",adjU[j][k]+1);
      printf("\n\n"); fflush(stdout);
#endif

  } // end for j


  
  free(idxLj);
  free(idxUj);
  free(pathlen); 
  free(flag);
  free(flaglist);
  free(queue);
  free(factorn);
  free(factord);

  
  
  startblock=endblock=0; // initial beginning of the current block
  cnt_fill=0; // counter for additional fill-in when merging blocks
  *nblocks=0; // number of diagonal blocks of A
  cntL=cntU=0;   // counter for accumulated index lists
  for (j=0; j<n; j++) {
      
#ifdef PRINT_INFO1
      printf("startblock=%6ld, endblock=%6ld\n",startblock,endblock);
      for (m=0; m<cntL; m++)
	  printf("%6ld",idxL[m]);
      printf("\n");
      for (m=0; m<cntU; m++)
	  printf("%6ld",idxU[m]);
      printf("\n");
      printf("current column/row j=%6ld\n",j);
      fflush(stdout);
#endif

      idxLj=adjL[j];
      cntLj=nzL[j];
      idxUj=adjU[j];
      cntUj=nzU[j];
      
      // we just started this block
      if (j==startblock) {
#ifdef PRINT_INFO1
	 printf("start block with column %ld\n",j);
#endif
	 // init idxL, idxU
	 cntL=cntLj;
	 memcpy(idxL,idxLj,(size_t)cntLj*sizeof(integer));
	 cntU=cntUj;
	 memcpy(idxU,idxUj,(size_t)cntUj*sizeof(integer));
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
	 // swap merged list with idxL
	 cntL=l;
	 pi=idxL; idxL=buff; buff=pi;

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
	 // swap merged list with idxU
	 cntU=l;
	 pi=idxU; idxU=buff; buff=pi;

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

	    // init idxL, idxU
	    cntL=cntLj;
	    memcpy(idxL,idxLj,(size_t)cntLj*sizeof(integer));
	    cntU=cntUj;
	    memcpy(idxU,idxUj,(size_t)cntUj*sizeof(integer));
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
  
#ifdef PRINT_INFO
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
  printf("number of blocks: %ld, maximum block size: %ld, average %4.1lf(+-%8.1le)\n",*nblocks,r,rval,weightd);
#endif


  // release memory
  free(idxL);
  free(idxU);

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
      if (adjU[i]!=NULL)
	 free(adjU[i]);
  } // end for i
  free(adjL);
  free(adjU);
  free(nzL);
  free(nzU);

  
  return (n-j);
} // end ilupt_blocks

 

 
