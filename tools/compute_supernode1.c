/*  Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        $Id: compute_supernode.c 4146 2018-04-19 11:55:16Z bolle $ 

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



// #define PRINT_INFO
// #define PRINT_INFO1
// #define printf mexPrintf

// #define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

void compute_union(integer **ka, integer *nka,
		   integer *ja,  integer nja,
		   integer pos, integer *idxpos, integer flag_sorted);


integer SUPERNODES(SPARSEMATRIX *A,
		   REALS *SL, REALS *SR, integer *p, integer *invq,
		   integer *blocksize, integer *nblocks,
		   REALS droptol)
{  
   integer n=A->nc, *idx, *idxpos, cnt,i,j,k,l,m, ii,jj,*ia, *ja,*ka, *parent,
           **reachable, *nreachable, flag, flagI,par, startpos, *idxpos2;
   SPARSEMATRIX BL,BU;
#ifdef _PROFILING_
   double timeBegin,
          time_total=0.0,
          time_graph=0.0,
          time_etree=0.0,
          time_super=0.0;

   time_total=omp_get_wtime();
#endif

   // auxiliary array for indices and its inverse map
   idx    =(integer *)malloc((size_t)n*sizeof(integer));
   idxpos =(integer *)calloc((size_t)n,sizeof(integer));
   idxpos2=(integer *)calloc((size_t)n,sizeof(integer));
   // parent structure in the elimination tree of |A(q,p)|+|A(q,p)|^T
   parent=(integer *)calloc((size_t)n,sizeof(integer));
   // list of reachable nodes for each node and their number
   reachable =(integer **)malloc((size_t)n*sizeof(integer *));
   nreachable=(integer *) calloc((size_t)n,sizeof(integer));
   for (j=0; j<n; j++)
       reachable[j]=NULL;

#ifdef _PROFILING_
  timeBegin = omp_get_wtime();
#endif
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
   // use BL.ncol as counter
   ia=BL.ncol;
   for (j=0; j<n; j++) {
       l=A->ncol[p[j]];
       ja=A->rowind[p[j]];
       for (k=0; k<l; k++) {
	   i=invq[ja[k]];
	   // increment number of nonzeros in column i of B
	   ia[i]++;
       } // end for k
   } // end for j
   for (j=0; j<n; j++) {
       l=A->ncol[p[j]]+ia[j];
       BL.rowind[j]=(integer *)malloc((size_t)l*sizeof(integer));
       // reset counter
       ia[j]=0;
   } // end for j
   // copy A^T
   for (j=0; j<n; j++) {
       l=A->ncol[p[j]];
       ja=A->rowind[p[j]];
       for (k=0; k<l; k++) {
	   i=invq[ja[k]];
	   // supplement BL with pattern of A^T
	   m=ia[i];
	   BL.rowind[i][m]=j;
	   // increment number of nonzeros in column i of B
	   ia[i]=m+1;
       } // end for k
   } // end for j
   // merge A+A^T
   for (j=0; j<n; j++) {
       // check-mark entries of B=A^T
       l=BL.ncol[j];
       ja=BL.rowind[j];
       for (k=0; k<l; k++) {
	   i=ja[k];
	   idxpos[i]=k+1;
       } // end for k
       // check how many additional entries are required for |A(q,p)|+|A(q,p)|^T
       m=l;
       l=A->ncol[p[j]];
       ka=A->rowind[p[j]];
       for (k=0; k<l; k++) {
	   i=invq[ka[k]];
	   // entry does not exist yet in B=A^T
	   if (!idxpos[i])
	      m++;
       } // end for k
       // increase memory for nonzeros in B(:,j)
       BL.rowind[j]=(integer *)realloc(ja,(size_t)m*sizeof(integer));
       ja=BL.rowind[j];
       // now insert entries
       m=BL.ncol[j];
       for (k=0; k<l; k++) {
	   i=invq[ka[k]];
	   // entry does not exist yet
	   if (!idxpos[i])
	      ja[m++]=i;
       } // end for k
       // clear check mark array
       l=BL.ncol[j];
       for (k=0; k<l; k++) {
	   i=ja[k];
	   idxpos[i]=0;
       } // end for k
       // update number of nonzeros in B
       BL.ncol[j]=m;

       // sort entries in increasing order and use idx as buff
       QQSORTI(ja,idx,&m);

       // where does the strict lower triangular part start?
       startpos=m;
       for (k=0; k<m; k++) {
	   i=ja[k];
	   if (i>j) {
	      startpos=k;
	      break;
	   } // end if
       } // end for k

       // strict lower triangular part used for supernodes
       // copy strict lower triangular part
       BL.ncol[j]=m-startpos;
       BL.rowind[j]=(integer *)malloc((size_t)(m-startpos)*sizeof(integer));
       memcpy(BL.rowind[j],ja+startpos,(m-startpos)*sizeof(integer));
       
       // upper triangular part used for ETREE0
       // shorten and re-allocate column
       BU.ncol[j]=startpos;
       BU.rowind[j]=(integer *)realloc(ja,(size_t)startpos*sizeof(integer));
   } // end for j
   // PRINTMATRIX(&BU);
   // PRINTMATRIX(&BL);
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
   printf("computed supernode sizes\n");
   for (l=0; l<*nblocks; l++)
       printf("%4ld",blocksize[l]);
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

#ifdef _PROFILING_
   time_total=omp_get_wtime()-time_total;
   printf("profiling summary\n");
   printf("1) computation of A+A^T                       %12.4le\n",time_graph);
   printf("2) computation elimination tree               %12.4le\n",time_etree);
   printf("3) detection of supernodal structure          %12.4le\n",time_super);
   printf("Total SUPERNODE time %12.4le\n\n",time_total);

   fflush(stdout);
#endif
   
   return (n-j);
} // end compute_supernode


#ifdef _DOUBLE_REAL_
// compute union of two index sets
// ka, nka destination index set and its initial size, to be increased
// a negative value for nka indicates that the array ka was simply shifted
// and does not necessarily fulfill that its entries are greater than pos
// ja, nja second index set and its size	    
// pos is a position filter to restrict ja to entries greater than pos
// idxpos is sufficiently large array of size at least nka, initialized
// with 0 prior to the call of this function and set back to 0 on return
void compute_union(integer **ka, integer *nka,
		   integer *ja,  integer nja,
		   integer pos, integer *idxpos, integer flag_sorted)
{
  integer *idx=*ka,*idx2,l, cnt, i,m, count_number_compressed=-1;

#ifdef PRINT_INFO
   printf("compute union, merge up to %4ld nodes >%4ld\n", *nka+nja, pos+1);
   fflush(stdout);
#endif
   // destination is void, simply copy/move entries
   if (*nka==0) {
      // is there anything to copy/move?
      if (nja>0) {
	 m=ja[0];
	 // array
	 if (flag_sorted && m>pos) {
	    // move pointer
#ifdef PRINT_INFO1
	    printf("move pointers\n");fflush(stdout);
#endif
	    *ka=ja; *nka=nja;
	    // do not free ja!
	 }
	 else {
	    // simply shift entries and move pointer
	    if (flag_sorted) {
	       m=nja;
	       for (l=0; l<nja; l++) {
		   i=ja[l];
		   if (i>pos) {
		      m=l;
		      break;
		   }
	       } // end for l
	       cnt=0;
	       for (l=m; l<nja; l++)
		   ja[cnt++]=ja[l];
	    }
	    else {
	       cnt=0;
	       for (l=0; l<nja; l++) {
	           i=ja[l];
		   if (i>pos)
		      ja[cnt++]=i;
	       } // end for l
	    }
	    // now move pointer and re-allocate
	    *ka=ja; *nka=cnt;
	    *ka=(integer *)realloc(ja,(size_t)cnt*sizeof(integer));
	 } // end if flag_sorted and m>pos
      } // end if nja>0
      else {
#ifdef PRINT_INFO1
	 printf("nothing to do\n");fflush(stdout);
#endif
      }
   } // end if *nka==0
   else if (nja==0) { // nothing to be added, keep values
      return;
   } // end if-elsif nja=0
   else { // nka!=0 and nja>0 a real merger is required
      // so far the left array *ka was only shifted from a single child
      if (*nka<0) {
	 *nka=-*nka;
	 // compress entries of *ka and check mark the remaining ones
	 idx2=idx;
	 m=0;
	 for (l=0; l<*nka; l++) {
	     i=idx[l];
	     if (i>pos) {
	        idx2[m++]=i;
	        idxpos[i]=m;
	     }
	 } // end for l
	 *nka=m;
	 // store how many of the shifted entries are left over
	 count_number_compressed=m;
      } // end if *nka<0
      else { // the array *ka is already properly filtered
	 // only check mark entries of *ka
	 for (l=0; l<*nka; l++) 
	     idxpos[idx[l]]=l+1;
      }
      // next check how many additional entries will enter
      if (flag_sorted) {
	 // detect position m where all indices >pos start
	 m=nja;
	 for (l=0; l<nja; l++) {
	     i=ja[l];
	     if (i>pos) {
	        m=l;
		break;
	     }
	 } // end for
	 cnt=0;
	 for (l=m; l<nja; l++) {
	     i=ja[l];
	     if (!idxpos[i])
	       cnt++;
	 } // end for l
      } // end if flag_sorted
      else {
	 cnt=0;
	 for (l=0; l<nja; l++) {
	     i=ja[l];
	     if (!idxpos[i] && i>pos)
	        cnt++;
	 } // end for
      } // end if-else flag_sorted
      if (cnt>0) {
	 *ka=(integer *)realloc(idx,(size_t)(*nka+cnt)*sizeof(integer));
	 // insert additional entries
	 idx=*ka+*nka;
	 if (flag_sorted) {
	    cnt=0;
	    for (l=m; l<nja; l++) {
	        i=ja[l];
		if (!idxpos[i])
		   idx[cnt++]=i;
	    } // end for
	 }
	 else {
	    cnt=0;
	    for (l=0; l<nja; l++) {
	        i=ja[l];
		if (!idxpos[i] && i>pos)
		   idx[cnt++]=i;
	    } // end for
	 } // end if-else flag_sorted
      } // end if cnt>0
      // uncheck old entries
      idx=*ka;
      for (l=0; l<*nka; l++) 
	idxpos[idx[l]]=0;
      *nka+=cnt;

      // nka used to be negative in advance, i.e. it refers to a shifted
      // set from a single child node
      if (count_number_compressed>=0) {
	 // Now check whether we simply compressed it or we augmented the
	 // set with new indices
	 if (*nka==count_number_compressed) {
	    printf("union is a simple compress to %4ld nodes\n", *nka);
	    fflush(stdout);
	    // *nka=-*nka;
	 } // end if
      } // end if
      
      // give up ja
      free(ja); ja=NULL; nja=0;
   } // end if-elseif-else *nka==0
#ifdef PRINT_INFO
   printf("union computed, %4ld nodes >%4ld merged\n", *nka, pos+1);
   fflush(stdout);
#endif
} // compute union
#endif
