/*  Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        $Id: spdapp_supernode.c 5123 2019-06-21 16:45:08Z bolle $ 

    Notice:

	Copyright (c) 2018 by TU Braunschweig.  All Rights Reserved.

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
#define MYABS(A)        (((A)>0)?(A):(-(A)))
#define STDERR          stderr
#define STDOUT          stdout


#define TWO_BY_TWO_THRESHOLD 0.1
#define TWO_BY_TWO_BOUND     1.5



// #define PRINT_INFO
// #define PRINT_INFO1
// #define PRINT_INFO2
// #define printf mexPrintf

// #define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

void compute_union(integer **ka, integer *nka,
		   integer *ja,  integer nja,
		   integer pos, integer *idxpos, integer flag_sorted);


integer SPDAPPSUPERNODES(SPARSEMATRIX *A,
			 REALS *SL, REALS *SR, integer *p, integer *invq,
			 integer *blocksize, integer *nblocks,
			 REALS droptol)
{  
   integer n=A->nc, *idx, *idxpos, cnt,i,j,k,l,m, ii,jj,kk, *ia, *ja,*ka, *parent,
           **reachable, *nreachable, flag, flagI,par, startpos, *idxpos2;
   SPARSEMATRIX BL,BU,BLe;
   FLOAT   val,          // temporary scalar numerical value
           *a;

   REALS   a11,a21,a22,   // auxiliary variables for determinant
           *Adiag,       // array of diagonal entries
           sdroptol=sqrt(droptol), srval,
           rval;         // temporary scalar real numerical value
   
#ifdef _PROFILING_
   double timeBegin,
          time_total=0.0,
          time_graph=0.0,
          time_etree=0.0,
          time_super=0.0;

   time_total=omp_get_wtime();
#endif


   
#ifdef _PROFILING_
   timeBegin = omp_get_wtime();
#endif
   // diagonal entries
   Adiag  =(REALS *)malloc((size_t)n*sizeof(REALS));

   for (i=0; i<n; i++)
       Adiag[i]=0.0;

   // extract diagonal entries from A and store them separately
   for (j=0; j<n; j++) {
       jj=p[j];
       l=A->ncol[jj];
       ja=A->rowind[jj];
       a=A->val[jj];
       for (k=0; k<l; k++) {
	   ii=ja[k];
	   i=invq[ii];
	   if (i==j) {
	      if (SL==NULL) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Adiag[i]=sqrt(fabs(a[k]));
#else
		 Adiag[i]=sqrt(fabs(a[k].r));
#endif
	      }
	      else {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Adiag[i]=sqrt(fabs(SL[jj]*a[k]*SL[jj]));
#else
		 Adiag[i]=sqrt(fabs(SL[jj]*a[k].r*SL[jj]));
#endif
	      }
	      break;
	   } // end if i=j
       } // end for k
   } // end for j
   // diagonal entries extracted, strictly speaking: their square root

   // now compute the undirected graph

   // auxiliary array for indices and its inverse map
   idx    =(integer *)calloc((size_t)n,sizeof(integer));
   idxpos =(integer *)calloc((size_t)n,sizeof(integer));
   idxpos2=(integer *)calloc((size_t)n,sizeof(integer));
   // parent structure in the elimination tree of |A(q,p)|+|A(q,p)|^T
   parent=(integer *)calloc((size_t)n,sizeof(integer));
   // list of reachable nodes for each node and their number
   reachable =(integer **)malloc((size_t)n*sizeof(integer *));
   nreachable=(integer *) calloc((size_t)n,sizeof(integer));
   for (j=0; j<n; j++)
       reachable[j]=NULL;

   // compute B as pattern of |A(q,p)|+|A(q,p)|^T and make sure that its
   // row entries are sorted in increasing order
   // starting position of the strict lower triangular part
   BL.nr=BL.nc=n;
   BL.ncol  =(integer *) calloc((size_t)n,sizeof(integer));
   BL.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
   BL.val   =NULL;
   BU.nr=BU.nc=n;
   BU.ncol  =(integer *) calloc((size_t)n,sizeof(integer));
   BU.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
   BU.val   =NULL;

   // Four-pass algorithm to compute the graph of
   // |A(p,p)|+|A(p,p)|^T 
   // first pass:  determine memory requirement for graph of |A(p,p)|^T
   // second pass: determine memory requirement for graph of |A(p,p)|
   //              and allocate memory
   // third pass:  copy indices of the graph of |A(p,p)|^T
   // fourth pass: copy remaining indices of the graph of |A(p,p)|,
   //              sort indices and split graph into its upper
   //              triangular part (used for the elimination tree)
   //              as well as its strict lower triangular part
   //              (needed by supernodal detection)
   
   // first pass:  determine memory requirement for block graph of |A(p,p)|^T
   // use BL.ncol as counter
   ia=BL.ncol;
   // use idx to count the large entries in each column of A
   for (j=0; j<n; j++) {
       // short cuts column p[j]
       l=A->ncol[p[j]];
       ja=A->rowind[p[j]];
       a=A->val[p[j]];
       rval =droptol *Adiag[j];
       srval=sdroptol*Adiag[j];
       // for all indices i of A(p,p(j)) increment associated row
       // counter in the graph
       for (k=0; k<l; k++) {
	   // row index i in A(p(i),p(j))
	   i=invq[ja[k]];
	   // increment number of nonzeros in column i of B
	   // filter small entries a priori
	   if (FABS(a[k])>=rval*Adiag[i]) {
	      ia[i]++;
	      // also increase number of large nonzeros in A(p,p(j))
	      idx[j]++;
	   } // end if
       } // end for k
   } // end for j

   // second pass: determine memory requirement for graph of |A(p,p)|+A(p,p)|^T
   //              and allocate memory
   for (j=0; j<n; j++) {
       // total number of large nonzeros in column j
       l=idx[j]+ia[j];
       // allocate sufficient memory for block column j
       BL.rowind[j]=(integer *)malloc((size_t)l*sizeof(integer));
       // reset counter
       ia[j]=0;
   } // end for j
   // third pass:  copy indices of the block graph of |A(p,p)|^T
   for (j=0; j<n; j++) {
       // short cuts column p[j]
       l=A->ncol[p[j]];
       ja=A->rowind[p[j]];
       a=A->val[p[j]];
       rval=droptol*Adiag[j];
       srval=sdroptol*Adiag[j];
       // transfer indices directly
       for (k=0; k<l; k++) {
	   // row index i in A(p(i),p(j))
	   i=invq[ja[k]];
	   // filter small entries a priori
	   if (FABS(a[k])>=rval*Adiag[i]) {
	      // counter for nonzeros in column i of B
	      m=ia[i];
	      // supplement BL with pattern of A^T
	      if (FABS(a[k])>=srval*Adiag[i] || j<=i)
		 BL.rowind[i][m]=j;
	      else
		 // flag smaller, strict lower triangular indices negative but
		 // keep in mind that indices start with 0, thus we even need
		 // to shift by -1
		 BL.rowind[i][m]=-j-1;
	      // BL.rowind[i][m]=j;
	      // increment number of nonzeros in column i of B
	      ia[i]=m+1;
	   } // end if 
       } // end for k
   } // end for j

   
   // fourth pass: copy remaining indices of the graph of |A(p,p)|
   //              we thus merge the graphs of A(p,p) and A(p,p)^T
   //              also sort indices and split graph into its upper
   //              triangular part (used for the elimination tree)
   //              as well as its strict lower triangular part
   //              (needed by supernodal detection)
   for (j=0; j<n; j++) {
       // transfer remaining entries from A(p,p(j)) directly
       // check mark existing entries in B(:,j)
       // short cuts column j of B
       l=BL.ncol[j];
       ja=BL.rowind[j];
       for (k=0; k<l; k++) {
	   i=ja[k];
	   // remember that smaller indices are negative and shifted
	   if (i<0)
	      i=-i-1;
	   idxpos[i]=k+1;
       } // end for k
       // check how many additional entries are required for |A(q,p)|+|A(q,p)|^T
       // current number of nonzeros
       m=l;
       // short cuts column p[j]
       l=A->ncol[p[j]];
       ka=A->rowind[p[j]];
       a=A->val[p[j]];
       rval=droptol*Adiag[j];
       for (k=0; k<l; k++) {
	   // row index i of A(p(i),p(j))
	   i=invq[ka[k]];
	   // entry does not exist yet in B=A^T
	   // also filter small entries a priori
	   if (!idxpos[i] && FABS(a[k])>=rval*Adiag[i])
	      m++;
       } // end for k
       // increase memory for nonzeros in B(:,j)
       BL.rowind[j]=(integer *)realloc(ja,(size_t)m*sizeof(integer));
       // update short cut
       ja=BL.rowind[j];
       // now insert entries
       // number of existing entries
       m=BL.ncol[j];
       for (k=0; k<l; k++) {
	   // row index i of A(p(i),p(j))
	   i=invq[ka[k]];
	   // entry does not exist yet
	   // also filter small entries a priori
	   if (!idxpos[i] && FABS(a[k])>=rval*Adiag[i]) {
	      if (FABS(a[k])>=srval*Adiag[i] || i<=j)
		 ja[m++]=i;
	      else
		 // flag smaller indices negative but keep in mind that indices
		 // start with 0, thus we even need to shift by -1
	      ja[m++]=-i-1;
	      // ja[m++]=i;
	   }
       } // end for k
       // clear check mark array
       l=BL.ncol[j];
       for (k=0; k<l; k++) {
	   // row index i in the graph
	   i=ja[k];
	   // remember that smaller indices are negative and shifted
	   if (i<0)
	      i=-i-1;
	   // clear flag
	   idxpos[i]=0;
       } // end for k
       // update number of nonzeros in B
       BL.ncol[j]=m;

       // now sort entries and split B into upper and strict lower
       // triangular part

       // sort entries in increasing order and use idx as buff
       QQSORTI(ja,idx,&m);


       // where does the upper triangular part start?
       startpos=m;
       for (k=0; k<m; k++) {
	   i=ja[k];
	   if (i>=0) {
	      startpos=k;
	      break;
	   } // end if
       } // end for k
       // where does the strict lower triangular part start?
       l=m;
       for (; k<m; k++) {
	   i=ja[k];
	   if (i>j) {
	      l=k;
	      break;
	   } // end if
       } // end for k
       // flip sign back and shift
       for (k=0; k<startpos; k++) {
	   i=ja[k];
	   ja[k]=-i-1;
       } // end for k

       
       // strict lower triangular part with larger entries used for supernodes
       // copy strict lower triangular part with larger entries
       BL.ncol[j]=m-l;
       BL.rowind[j]=(integer *)malloc((size_t)(m-l)*sizeof(integer));
       memcpy(BL.rowind[j],ja+l,(m-l)*sizeof(integer));

       // upper triangular part used for ETREE0
       // copy strict lower triangular part
       BU.ncol[j]=l-startpos;
       BU.rowind[j]=(integer *)malloc((size_t)(l-startpos)*sizeof(integer));
       memcpy(BL.rowind[j],ja+startpos,(l-startpos)*sizeof(integer));
       
       // strict lower triangular part with smaller indices also used for supernodes
       // shorten and re-allocate column
       BLe.ncol[j]=startpos;
       BLe.rowind[j]=(integer *)realloc(ja,(size_t)startpos*sizeof(integer));
   } // end for j

#ifdef PRINT_INFO2
   i=0;
   for (j=0; j<n; j++)
       i+=BU.ncol[j];
   printf("BU, number of nonzeros: %ld\n",i);
   i=0;
   for (j=0; j<n; j++)
       i+=BL.ncol[j];
   printf("BL, number of nonzeros: %ld\n",i);
   fflush(stdout);
   i=0;
   for (j=0; j<n; j++)
       i+=BLe.ncol[j];
   printf("BLe, number of nonzeros: %ld\n",i);
   fflush(stdout);
#endif
   // PRINTMATRIX(&BU);
   // PRINTMATRIX(&BL);
   // PRINTMATRIX(&BLe);
   // end computation of B as pattern of |A(q,p)|+|A(q,p)|^T
#ifdef _PROFILING_
   time_graph=omp_get_wtime()-timeBegin;
#endif

   
   // compute elimination tree of |A(q,p)|+|A(q,p)|^T
#ifdef _PROFILING_
   timeBegin = omp_get_wtime();
#endif
   ETREE0(&BU,parent,idx);
   SPARSEDELETE(&BU);
#ifdef _PROFILING_
   time_etree=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_INFO
   printf("elimination tree, parent array\n"); fflush(stdout);
   for (j=0; j<n; j++)
       printf("%4ld",parent[j]+1);
   printf("\n"); fflush(stdout);
#endif

   
#ifdef _PROFILING_
   timeBegin = omp_get_wtime();
#endif
   // compute supernodal structure and loop over all nodes
   cnt=0;
   j=0;
   *nblocks=-1;
   flagI=-1;
   while (j<n) {
         // start a new supernode beginning with node j
         (*nblocks)++;
         blocksize[*nblocks]=1;
#ifdef PRINT_INFO1
	 printf("new supernode %4ld, start with node %ld\n", *nblocks+1, j+1);
	 fflush(stdout);
#endif
	 
	 // nodes that are reachable from j
	 // check whether we resume from a previous step
	 // NO!
	 if (flagI) {
	    // 1. new nodes reachable via G(A), indices i>j
	    ja=BL.rowind[j];
	    m=BL.ncol[j];
	    // 2. parent node from the elimination tree
	    par=parent[j];
	    // build union with existing nodes >j inherited from its descendants
	    // and give up BL.rowind[j]
#ifdef PRINT_INFO1
	    printf("parent[%ld]=%ld\n",j+1,par+1);
	    if (nreachable[j]<0)
	       printf("compress and unite reachable[%ld] with its neighbours\n",j+1);
	    else
	       printf("unite reachable[%ld] with its neighbours\n",j+1);
	    fflush(stdout);
#endif
	    compute_union(reachable+j,nreachable+j, ja,m, j, idxpos2, 1);
	    BL.rowind[j]=NULL; BL.ncol[j]=0;
	 } // end if 
	 else
	    flagI=-1;
#ifdef PRINT_INFO1
	 if (nreachable[j]<0)
	    printf("(shifted) reachable set node %ld\n",j+1);
	 else
	    printf("reachable set node %ld\n",j+1);
	 for (l=0; l<MYABS(nreachable[j]); l++)
	     printf("%4ld",reachable[j][l]+1);
	 printf("\n");
	 fflush(stdout);
#endif

	 // can we continue to build up a supernode?
	 // yes
	 if (par==j+1) {
	    // init supernode with reachable set of j
	    cnt=MYABS(nreachable[j]);
	    memcpy(idx,reachable[j],cnt*sizeof(integer));
	    // now check mark nodes of the reachable set of j
	    for (l=0; l<cnt; l++) 
		idxpos[idx[l]]=l+1;
	    // always move reachable nodes to its parent
	    // to merge these arrays, first check mark existing entries
	    // and give up old data
	    if (nreachable[par]==0) {
#ifdef PRINT_INFO1
	       printf("shift reachable set from %ld to %ld\n",j+1,par+1); fflush(stdout);
#endif
	       // only shift and indicate shift via negative value
	       reachable[par] = reachable[j];
	       nreachable[par]=-MYABS(nreachable[j]);
	    }
	    else {
#ifdef PRINT_INFO1
	       if (nreachable[par]<0)
		  printf("compress and unite reachable[%ld] with reachable[%ld]\n",par+1,j+1);
	       else
		  printf("unite reachable[%ld] with reachable[%ld]\n",par+1,j+1);
	       fflush(stdout);
#endif
	       compute_union(reachable+par,nreachable+par,
			     reachable[j],MYABS(nreachable[j]), par, idxpos2,0);
#ifdef PRINT_INFO1
	       if (nreachable[par]<0)
		  printf("(shifted) reachable set node %ld\n",par+1);
	       else
		  printf("reachable set node %ld\n",par+1);
	       for (l=0; l<MYABS(nreachable[par]); l++)
		   printf("%4ld",reachable[par][l]+1);
	       printf("\n");
	       fflush(stdout);
#endif
	    }
	    reachable[j]=NULL; nreachable[j]=0;
	    

	    // move along the elimination tree as long as the parent nodes
	    // have consecutive numbers 
	    flag=-1;
	    k=j+1;
	    while (flag) {
	          // 1. new nodes
	          ja=BL.rowind[k];
		  m=BL.ncol[k];
		  //  2. parent
		  par=parent[k];
		  // build union with existing nodes >k inherited from its
		  // descendants
		  // and give up BL.rowind[k]
		  if (nreachable[k]<0) {
		     // only check whether indices of BL.rowind[k] are already
		     // in use 
		     ii=0;
		     for (l=0; l<m; l++)
		         if (idxpos[ja[l]]==0) {
			    ii=-1;
			    break;
			 } // end if
		     // NO, they are not!
		     if (ii) {
#ifdef PRINT_INFO1
		        printf("compress and unite reachable[%ld] with its neighbours\n",k+1);
			fflush(stdout);
#endif
		        // compress and unite
		        compute_union(reachable+k,nreachable+k, ja,m, k, idxpos2,1);
		     }
		     else 
		        free(ja);
		     ja=reachable[k];
		     m=nreachable[k];
		  }
		  else { // nreachable[k]>=0
		     // unite and check all
#ifdef PRINT_INFO1
		     printf("unite reachable[%ld] with its neighbours\n",k+1);
		     fflush(stdout);
#endif
		     compute_union(reachable+k,nreachable+k, ja,m, k, idxpos2,1);
		     // check whether new reachable indices are already in use
		     ja=reachable[k];
		     m=nreachable[k];
		     ii=0;
		     for (l=0; l<m; l++)
		         if (idxpos[ja[l]]==0) {
			    ii=-1;
			    break;
			 } // end if
		  }
		  BL.rowind[k]=NULL; BL.ncol[k]=0;
#ifdef PRINT_INFO1
		  if (m<0)
		     printf("(shifted) reachable set node %ld\n",k+1);
		  else
		     printf("reachable set node %ld\n",k+1);
		  for (l=0; l<MYABS(m); l++)
		      printf("%4ld",ja[l]+1);
		  printf("\n");
		  fflush(stdout);
#endif
		  
		  // YES, they are!
		  if (ii==0) {
		     // no additional existing nodes >k
		     // add k to supernode nblocks
		     blocksize[*nblocks]++;
#ifdef PRINT_INFO1
		     printf("supernode %4ld, add node %ld, current size=%4ld\n", *nblocks+1, k+1,blocksize[*nblocks]);
		     fflush(stdout);
#endif
		     // always move reachable nodes to its parent
		     if (par>=0) {
		        // ... and give up old data
		        if (nreachable[par]==0) {
#ifdef PRINT_INFO1
			   printf("shift reachable set from %ld to %ld\n",k+1,par+1); fflush(stdout);
#endif
			   // only shift and indicate shift via negative value
			   reachable[par] = ja;
			   nreachable[par]=-MYABS(m);
			}
			else {
#ifdef PRINT_INFO1
			   if (nreachable[par]<0)
			      printf("compress and unite reachable[%ld] with reachable[%ld]\n",par+1,k+1);
			   else
			      printf("unite reachable[%ld] with reachable[%ld]\n",par+1,k+1);
			   fflush(stdout);
#endif
			   compute_union(reachable+par,nreachable+par,
					 ja,MYABS(m), par, idxpos2,0);
#ifdef PRINT_INFO1
			   if (nreachable[par]<0)
			      printf("(shifted) reachable set node %ld\n",par+1);
			   else
			      printf("reachable set node %ld\n",par+1);
			   for (l=0; l<MYABS(nreachable[par]); l++)
			       printf("%4ld",reachable[par][l]+1);
			   printf("\n");
			   fflush(stdout);
#endif
			}
		     }
		     else {
		        // only give up old data
		        free(reachable[k]);
		     }
		     reachable[k]=NULL; nreachable[k]=0;

		     // is there another potential supernode available?
		     // yes
		     if (par==k+1)
		        k=k+1;
		     else { // no, k is the last one in this sequence
		        flag=0;
			for (l=0; l<cnt; l++) 
			    idxpos[idx[l]]=0;
			cnt=0;
			j=k+1;
		     } // end if-else
		  } // end if ii=0
		  else { // NO, a new supernode is going to start with j=k
#ifdef PRINT_INFO1
		     printf("terminate supernode %4ld, and restart next supernode with node %ld\n", *nblocks+1, k+1);
		     fflush(stdout);
#endif
		     flag=0;
		     for (l=0; l<cnt; l++) 
		         idxpos[idx[l]]=0;
		     cnt=0;
		     j=k;
		     flagI=0;
		  } // end if-else
	    } // end while flag
	 } // end if par=j+1
	 else { // NO,  the supernode remains scalar
#ifdef PRINT_INFO1
	    printf("supernode %4ld remains scalar with node %ld\n", *nblocks+1, j+1);
	    fflush(stdout);
#endif
	   // always move reachable nodes to its parent
	   if (par>=0) {
	      // ... and give up old data
	      if (nreachable[par]==0) {
		 // only shift and indicate shift via negative value
#ifdef PRINT_INFO1
		 printf("shift reachable set from %ld to %ld\n",j+1,par+1); fflush(stdout);
#endif
		 reachable[par] = reachable[j];
		 nreachable[par]=-MYABS(nreachable[j]);
	      }
	      else {
#ifdef PRINT_INFO1
		 if (nreachable[par]<0)
		    printf("compress and unite reachable[%ld] with reachable[%ld]\n",par+1,j+1);
		 else
		    printf("unite reachable[%ld] with reachable[%ld]\n",par+1,j+1);
		 fflush(stdout);
#endif
		 compute_union(reachable+par,nreachable+par,
			       reachable[j],MYABS(nreachable[j]), par, idxpos2,0);
#ifdef PRINT_INFO1
		 if (nreachable[par]<0)
		    printf("(shifted) reachable set node %ld\n",par+1);
		 else
		    printf("reachable set node %ld\n",par+1);
		 for (l=0; l<MYABS(nreachable[par]); l++)
		     printf("%4ld",reachable[par][l]+1);
		 printf("\n");
		 fflush(stdout);
#endif
	      }
	   }
	   else {
	      // only give up old data
	      free(reachable[j]); 
	   }
	   reachable[j]=NULL; nreachable[j]=0;

	   // advance to next column
	   j=j+1;
	 } // end if-else
   } // end while j<n
   // finally adjust the number of blocks
   (*nblocks)++;
#ifdef _PROFILING_
   time_super=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_INFO
   printf("computed supernode sizes\n");
   for (l=0; l<*nblocks; l++)
       printf("%6ld",blocksize[l]);
   printf("\n"); fflush(stdout);
#endif
   
   SPARSEDELETE(&BL);

   
   free(idx);
   free(idxpos);
   free(idxpos2);
   free(parent);
   for (l=0; l<n; l++)
       if (nreachable[l]!=0)
	  free(reachable[l]);
   free(reachable);
   free(nreachable);
   free(Adiag);
   
#ifdef _PROFILING_
   time_total=omp_get_wtime()-time_total;
   printf("profiling summary\n");
   printf("1) computation of A+A^T                       %12.4le\n",time_graph);
   printf("2) computation elimination tree               %12.4le\n",time_etree);
   printf("3) detection of supernodal structure          %12.4le\n",time_super);
   printf("Total SUPERNODE time %12.4le\n\n",time_total);

   fflush(stdout);
#endif
   
   return (0);
} // end MYSYMSUPERNODES


