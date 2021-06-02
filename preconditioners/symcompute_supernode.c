/*  Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        $Id: symcompute_supernode.c 4087 2018-03-19 18:58:22Z bolle $ 

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


#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYSYMSUPERNODES     SSYMsupernodes
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYSYMSUPERNODES     DSYMsupernodes
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYSYMSUPERNODES     CSYMsupernodes
#define CONJG(A)       (A)
#else // double complex
#define MYSYMSUPERNODES     ZSYMsupernodes
#define CONJG(A)       (A)
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYSYMSUPERNODES     SSYMsupernodes
#define CONJG(A)       (A)
#elif defined _DOUBLE_REAL_
#define MYSYMSUPERNODES     DSYMsupernodes
#define CONJG(A)       (A)
#elif defined _SINGLE_COMPLEX_
#define MYSYMSUPERNODES     CHERsupernodes
#define CONJG(A)       (-(A))
#else // double complex
#define MYSYMSUPERNODES     ZHERsupernodes
#define CONJG(A)       (-(A))
#endif //-elif-else single-real

#endif //-else complex-symmetric



#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define MYABS(A)        (((A)>0)?(A):(-(A)))
#define STDERR          stderr
#define STDOUT          stdout


#define TWO_BY_TWO_THRESHOLD 0.1
#define TWO_BY_TWO_BOUND     1.5



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


integer MYSYMSUPERNODES(SPARSEMATRIX *A,
			REALS *SL, REALS *SR, integer *p, integer *invq,
			integer *blocksize, integer *nblocks,
			REALS droptol)
{  
   integer n=A->nc, *idx, *idxpos, cnt,i,j,k,l,m, ii,jj,kk, *ia, *ja,*ka, *parent,
           *blockstructure, colj,
           **reachable, *nreachable, flag, flagI,par, startpos, *idxpos2;
   SPARSEMATRIX BL,BU;
   FLOAT   val,          // temporary scalar numerical value
           *a,
          *Adiag,       // array of diagonal entries 
          *Asdiag;      // array of sub-diagonal entries

   REALS   a11,a21,a22,   // auxiliary variables for determinant
           bnd,bndl,bndr,
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
   Adiag  =(FLOAT *)malloc(n*sizeof(FLOAT));
   // sub-diagonal entries 
   Asdiag =(FLOAT *)malloc(n*sizeof(FLOAT));

   // block structure for 1x1 and 2x2 pivots
   blockstructure=(integer *)calloc(n+1,sizeof(integer));

#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
   for (i=0; i<n; i++)
       Adiag[i]=Asdiag[i]=0.0;
#else
   for (i=0; i<n; i++)
       Adiag[i].r=Asdiag[i].r=Adiag[i].i=Asdiag[i].i=0.0;
#endif

   // extract diagonal and sub-diagonal entries from A and store them separately
   for (j=0; j<n; j++) {
       jj=p[j];
       l=A->ncol[jj];
       ja=A->rowind[jj];
       a=A->val[jj];
       for (k=0; k<l; k++) {
	   ii=ja[k];
	   i=invq[ii];
	   if (i==j)
	      if (SL==NULL) 
		 Adiag[i]=a[k];
	      else {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Adiag[i]=SL[jj]*a[k]*SL[jj];
#else
		 Adiag[i].r=SL[jj]*a[k].r*SL[jj];
		 Adiag[i].i=SL[jj]*a[k].i*SL[jj];
#endif
	      }
	   if (i==j+1)
	      if (SL==NULL) 
		 Asdiag[i]=a[k];
	      else {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Adiag[i]=SL[ii]*a[k]*SL[jj];
#else
		 Adiag[i].r=SL[ii]*a[k].r*SL[jj];
		 Adiag[i].i=SL[ii]*a[k].i*SL[jj];
#endif
	      }
	   if (i==j-1)
	      if (SL==NULL) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Asdiag[i]=a[k];
#else
		 Asdiag[i].r=a[k].r;
		 Asdiag[i].i=CONJG(a[k].i);
#endif
	      } // end if SL=0
	      else {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		 Asdiag[i]=SL[ii]*a[k]*SL[jj];
#else
		 Asdiag[i].r=SL[ii]*a[k].r*SL[jj];
		 Asdiag[i].i=SL[ii]*CONJG(a[k].i)*SL[jj];
#endif
	      }	// end if-else SL=0 
       } // end for k
   } // end for j
   // diagonal and sub-diagonal entries extracted

   // now determine when to prefer a 1x1 pivot and when to stick with a 2x2 pivot
   j=0;
   m=0;
   while (j<n) {

         // scalar bound
         bnd=FABS(Adiag[j]);
         // so far no additional column (i.e. no block approach)
	 colj=-2;

	 // does A(p[j-1],p[j-1]) already belong to a 2x2 pivot A(p[j-2:j-1],p[j-2:j-1])
	 flag=-1;
	 if (j>1)
	    // yes!
	    if (blockstructure[m]-blockstructure[m-1]==2)
	       flag=0;
   
	 // check whether a 1x1 pivot or a 2x2 pivot should be chosen
	 // first check preceding column j-1
         if (j>0 && flag) {
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
		  if (j>1) {
		     // yes, rescind previous 1x1 pivot decision
		     if (blockstructure[m]-blockstructure[m-1]==1) {
#ifdef PRINT_INFO1
		        printf("rescind 1x1 pivot %4ld, new 2x2 pivot (%4ld,%4ld)\n",j-1,j-1,j); fflush(stdout);
#endif
		        blockstructure[m]=blockstructure[m-1]+2;
		     }
		     else { // no, new 2x2 pivot
#ifdef PRINT_INFO1
		        printf("new 2x2 pivot (%4ld,%4ld)\n",j-1,j); fflush(stdout);
#endif
		        blockstructure[m+1]=blockstructure[m]+2;
			m++;
		     }
		  }
		  else { // new 2x2 pivot
#ifdef PRINT_INFO1
		     printf("new 2x2 pivot (%4ld,%4ld)\n",j-1,j); fflush(stdout);
#endif
		     blockstructure[m+1]=blockstructure[m]+2;
		     m++;
		  }
		  j++;
	       } // end if
	    } // end if a21 large enough
	 } // end if j>0

	 // possibly second check the subsequent column j+1
	 if (j<n-1 && colj<0) {
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
		  // new 2x2 pivot
#ifdef PRINT_INFO1
		  printf("new 2x2 pivot (%4ld,%4ld)\n",j,j+1); fflush(stdout);
#endif
		  blockstructure[m+1]=blockstructure[m]+2;
		  m++;
		  j+=2;
	       } // end if
	    } // end if a21 large enough
	 } // end if j<n-1

	 if (colj<0) {
	    // 1x1 pivot, possibly only temporarily
#ifdef PRINT_INFO1
	    printf("1x1 pivot %4ld, maybe temporarily\n",j); fflush(stdout);
#endif
	    blockstructure[m+1]=blockstructure[m]+1;
	    m++;
	    j++;
	 }
   } // end while j
#ifdef PRINT_INFO
   printf("in total we have %ld blocks\n",m);
   for (j=0; j<=m; j++) {
      printf("%6ld",blockstructure[j]);
   }
   printf("\n");
   fflush(stdout);
#endif

   free(Adiag);
   free(Asdiag);
   // computation of 1x1 pivots and 2x2 pivots complete, now compute
   // the undirected AND COMPRESSED graph

   
   
   // auxiliary array for indices and its inverse map
   idx    =(integer *)malloc((size_t)(n+1)*sizeof(integer));
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
   BL.nr=BL.nc=m;
   BL.ncol  =(integer *) calloc((size_t)m,sizeof(integer));
   BL.rowind=(integer **)malloc((size_t)m*sizeof(integer *));
   BL.val   =NULL;
   BU.nr=BU.nc=m;
   BU.ncol  =(integer *) calloc((size_t)m,sizeof(integer));
   BU.rowind=(integer **)malloc((size_t)m*sizeof(integer *));
   BU.val   =NULL;

   // compute inverse mapping w.r.t. blockstructure and use idx to represent
   // the inverse mapping
   for (i=0; i<m; i++) {
       // starting index of block i in the original matrix
       j=blockstructure[i];
       // index j belongs to block i
       idx[j]=i;
       // 2x2 block
       if (blockstructure[i+1]-blockstructure[i]==2)
	  // index j+1 also belongs to block i
          idx[j+1]=i;
   } // end for i


   // Five-pass algorithm to compute the compressed block graph of
   // |A(p,p)|+|A(p,p)|^T w.r.t. blockstructure
   // first pass:  determine memory requirement for block graph of |A(p,p)|^T
   // second pass: determine memory requirement for block graph of |A(p,p)|
   //              and allocate memory
   // third pass:  copy indices of the block graph of |A(p,p)|^T
   // fourth pass: copy remaining indices of the block graph of |A(p,p)|
   // fifth pass:  sort indices and split block graph into its upper
   //              triangular part (used for the elimination tree) and
   //              its strict lower triangular part (needed by supernodal
   //              detection)


   // first pass:  determine memory requirement for block graph of |A(p,p)|^T
   // use BL.ncol as counter
   ia=BL.ncol;
   // counter for the columns  of A
   j=0;
   // counter for the number of blocks
   m=0;
   while (j<n) {
         // short cuts column p[j]
         l=A->ncol[p[j]];
	 ja=A->rowind[p[j]];
	 // 2x2 pivot case, block column p(j) and p(j+1) contribute in common
	 if (blockstructure[m+1]-blockstructure[m]==2) {
	    // check mark row indices i of A(p,p(j)) and for any i increment
  	    // associated row counter in the block graph
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix B
		ii=idx[i];
		// flag entry
		idxpos[ii]=1;
		// increment number of nonzeros in column ii of B
		ia[ii]++;
	    } // end for k
	    // for all indices i of A(p,p(j+1)) increment associated row
	    // counter in the block graph if it has not been marked yet
	    // short cuts column p[j+1]
	    l=A->ncol[p[j+1]];
	    ja=A->rowind[p[j+1]];
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j+1))
	        i=invq[ja[k]];
		// row index in the compressed matrix
		ii=idx[i];
		// is the row index already checked by column p(j)?
		// NO!
		if (!idxpos[ii]) 
		   // increment number of nonzeros in column ii of B
		   ia[ii]++;
	    } // end for k
	    // clear check marks
	    // short cuts column p[j]
	    l=A->ncol[p[j]];
	    ja=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index in the compressed matrix
		ii=idx[i];
		// clear check mark
		idxpos[ii]=0;
	    } // end for k
	    // advance to the index after the next one of A(p,p)
	    j+=2;
	 }
	 else {
	    // for all indices i of A(p,p(j)) increment associated row
	    // counter in the block graph
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index in the compressed matrix
		ii=idx[i];
		// increment number of nonzeros in column ii of B
		ia[ii]++;
	    } // end for k
	    // advance to the next column index of A(p,p)
	    j++;
	 } // end if-else 2x2 case
	 // simultaneously advance to the next block column ass. with A(p,p)
	 m++;
   } // end while j<n
#ifdef PRINT_INFO1
   printf("1. pass completed\n");
#endif

   // second pass: determine memory requirement for block graph of |A(p,p)|
   //              and allocate memory
   // counter for the number of blocks
   m=0;
   // counter for the columns  of A
   j=0;
   while (j<n) {
	 // 2x2 pivot case, the nonzeros in block column m is the union of the
         // nonzeros in column j and j+1 of A(p,p)
	 if (blockstructure[m+1]-blockstructure[m]==2) {
	    // compute nnz of column p[j:j+1]
	    // short cuts column p[j]
	    l=A->ncol[p[j]];
	    ja=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// flag entry
		idxpos[ii]=1;
	    } // end for k
	    // current number of nonzeros by column p(j)
	    kk=l;
	    // additional fill by column p[j+1]
	    // short cuts column p[j+1]
	    l=A->ncol[p[j+1]];
	    ja=A->rowind[p[j+1]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// Is this an additional entry compared with column p(j)?
		// YES!
		if (!idxpos[ii]) 
		   // increment number of nonzeros in block column p[j:j+1] of A
		   kk++;
	    } // end for k
	    // clear check mark array
	    // short cuts column p[j]
	    l=A->ncol[p[j]];
	    ja=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// clear check mark
		idxpos[ii]=0;
	    } // end for k

	    // total number of nonzeros in block column m
	    l=kk+ia[m];
	    // advance to the index after the next one of A(p,p)
	    j+=2;
	 }
	 else {
	    // total number of nonzeros in block column m
	    l=A->ncol[p[j]]+ia[m];
	    // advance to the next column index of A(p,p)
	    // advance to the next column index of A(p,p)
	    j++;
	 }

	 // allocate sufficient memory for block column m
	 BL.rowind[m]=(integer *)malloc((size_t)l*sizeof(integer));
	 // reset counter
	 ia[m]=0;
	 // simultaneously advance to the next block column ass. with A(p,p)
	 m++;
   } // end while j<n
#ifdef PRINT_INFO1
   printf("2. pass completed\n");
#endif

   // third pass:  copy indices of the block graph of |A(p,p)|^T
   // counter for the columns  of A
   j=0;
   // counter for the number of blocks
   m=0;
   while (j<n) {
         // short cuts column p[j]
         l=A->ncol[p[j]];
	 ja=A->rowind[p[j]];
	 // 2x2 pivot case, block column p(j) and p(j+1) contribute in common
	 if (blockstructure[m+1]-blockstructure[m]==2) {
	    // check mark row indices i of A(p,p(j)) and for any i transfer
  	    // associated row row index ii to the block graph
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index in the compressed matrix
		ii=idx[i];
		// flag entry
		idxpos[ii]=1;
		// counter for nonzeros in column ii of B
		jj=ia[ii]; 
		// supplement BL with pattern of A^T
		BL.rowind[ii][jj]=idx[j];
		// increment number of nonzeros in column ii of B
		ia[ii]=jj+1;
	    } // end for k
	    // for all indices i of A(p,p(j+1)) transfer associated row
	    // indices ii to the block graph if they have not been tranfered yet
	    // short cuts column p[j+1]
	    l=A->ncol[p[j+1]];
	    ja=A->rowind[p[j+1]];
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j+1))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// is the row index ii already inserted by column p(j)?
		// NO!
		if (!idxpos[ii]) {
		   // counter for nonzeros in column ii of B
		   jj=ia[ii]; 
		   // supplement BL with pattern of A^T
		   BL.rowind[ii][jj]=idx[j];
		   // increment number of nonzeros in column ii of B
		   ia[ii]=jj+1;
		}
	    } // end for k
	    // clear check marks
	    // short cuts column p[j]
	    l=A->ncol[p[j]];
	    ja=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index in the compressed matrix
		ii=idx[i];
		// clear check mark
		idxpos[ii]=0;
	    } // end for k
	    // advance to the index after the next one of A(p,p)
	    j+=2;
	 }
	 else {
	    // 1x1 pivot case, transfer ass. indices of the block graph directly
	    for (k=0; k<l; k++) {
	        // row index i in A(p(i),p(j))
	        i=invq[ja[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// counter for nonzeros in column ii of B
		jj=ia[ii]; 
		// supplement BL with pattern of A^T
		BL.rowind[ii][jj]=idx[j];
		// increment number of nonzeros in column ii of B
		ia[ii]=jj+1;
	    } // end for k
	    // advance to the next column index of A(p,p)
	    j++;
	 } // end if-else 2x2 case
	 // simultaneously advance to the next block column ass. with A(p,p)
	 m++;
   } // end while j<n
#ifdef PRINT_INFO1
   printf("3. pass completed\n");
#endif
   
   // fourth pass: copy remaining indices of the block graph of |A(p,p)|
   //              we thus merge the block graphs of A(p,p) and A(p,p)^T
   // counter for the columns  of A
   j=0;
   // counter for the number of blocks
   m=0;
   while (j<n) {
         // check-mark entries of B which refers to the block graph of A(p,p)^T
         // 2x2 pivot case
	 if (blockstructure[m+1]-blockstructure[m]==2) {
	    // check mark row ii indices of B(:,m)
	    // short cuts for column m of B
	    l=BL.ncol[m];
	    ja=BL.rowind[m];
	    for (k=0; k<l; k++) {
	        ii=ja[k];
		idxpos[ii]=k+1;
	    } // end for k
	    // check how many additional entries are required for |A(q,p)|+|A(q,p)|^T
	    // current number of nonzeros 
	    kk=l;
	    // short cuts for column j of A(p,p)
	    l=A->ncol[p[j]];
	    ka=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet in B
		if (!idxpos[ii]) {
		   // use -1 as check mark in order to distinguish between
		   // entries of B(:,m) and those induced by A(:,p(j))
		   idxpos[ii]=-1;
		   kk++;
		}
	    } // end for k
	    // short cuts for column j+1 of A(p,p)
	    l=A->ncol[p[j+1]];
	    ka=A->rowind[p[j+1]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j+1))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet neither in B(:,m) nor A(p,p(j))
		if (!idxpos[ii]) 
		   kk++;
	    } // end for k
	    // increase memory for nonzeros in B(:,j)
	    BL.rowind[m]=(integer *)realloc(ja,(size_t)kk*sizeof(integer));
	    // update short cut
	    ja=BL.rowind[m];
	    // now insert entries
	    // number of existing entries
	    kk=BL.ncol[m];
	    // short cuts column j of A(p,p)
	    l=A->ncol[p[j]];
	    ka=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet in B(:,m) but it could have been
		// flagged negative by a previous run through A(:,p(j))
		if (idxpos[ii]<=0)
		   ja[kk++]=ii;
	    } // end for k
	    // short cuts column j+1 of A(p,p)
	    l=A->ncol[p[j+1]];
	    ka=A->rowind[p[j+1]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j+1))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet, neither in B(:,m) nor in A(p,p(j))
		if (!idxpos[ii])
		   ja[kk++]=ii;
	    } // end for k
	    // clear check mark array
	    l=BL.ncol[m];
	    for (k=0; k<l; k++) {
		// row index ii in the compressed matrix
	        ii=ja[k];
		// clear flag
		idxpos[ii]=0;
	    } // end for k
	    l=A->ncol[p[j]];
	    ka=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// clear (additional) flag
		idxpos[ii]=0;
	    } // end for k
	    // advance to the index after the next one of A(p,p)
	    j+=2;
	 }
	 else {
	    // 1x1 pivot case transfer remaining entries from A(p,p(j)) directly
	    // check mark existing entries in B(:,m)
	    // short cuts column m of B
	    l=BL.ncol[m];
	    ja=BL.rowind[m];
	    for (k=0; k<l; k++) {
	        ii=ja[k];
		idxpos[ii]=k+1;
	    } // end for k
	    // check how many additional entries are required for |A(q,p)|+|A(q,p)|^T
	    // current number of nonzeros
	    kk=l;
	    // short cuts column p[j]
	    l=A->ncol[p[j]];
	    ka=A->rowind[p[j]];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet in B=A^T
		if (!idxpos[ii])
		   kk++;
	    } // end for k
	    // increase memory for nonzeros in B(:,j)
	    BL.rowind[m]=(integer *)realloc(ja,(size_t)kk*sizeof(integer));
	    // update short cut
	    ja=BL.rowind[m];
	    // now insert entries
	    // number of existing entries
	    kk=BL.ncol[m];
	    for (k=0; k<l; k++) {
	        // row index i of A(p(i),p(j))
	        i=invq[ka[k]];
		// row index ii in the compressed matrix
		ii=idx[i];
		// entry does not exist yet
		if (!idxpos[ii])
		   ja[kk++]=i;
	    } // end for k
	    // clear check mark array
	    l=BL.ncol[m];
	    for (k=0; k<l; k++) {
		// row index ii in the compressed matrix
	        ii=ja[k];
		// clear flag
		idxpos[ii]=0;
	    } // end for k
	    // advance to the next column index of A(p,p)
	    j++;
	 } // if-else 2x2 block case

	 // update number of nonzeros in B
	 BL.ncol[m]=kk;
	 // simultaneously advance to the next block column ass. with A(p,p)
	 m++;
   } // end while j<n
#ifdef PRINT_INFO1
   printf("4. pass completed\n");
#endif

   // fifth pass:  sort indices and split block graph into its upper
   //              triangular part (used for the elimination tree) and
   //              its strict lower triangular part (needed by supernodal
   //              detection)
   // split B into upper and strict lower triangular part
   for (jj=0; jj<m; jj++) {
       // short cuts column jj
       kk=BL.ncol[jj];
       ja=BL.rowind[jj];
       // sort entries in increasing order and use idx as buff
       QQSORTI(ja,idx,&kk);

       // where does the strict lower triangular part start?
       startpos=kk;
       for (l=0; l<kk; l++) {
	   ii=ja[l];
	   if (ii>jj) {
	      startpos=l;
	      break;
	   } // end if
       } // end for l
       // strict lower triangular part used for supernodes
       // copy strict lower triangular part
       BL.ncol[jj]=kk-startpos;
       BL.rowind[jj]=(integer *)malloc((size_t)(kk-startpos)*sizeof(integer));
       memcpy(BL.rowind[jj],ja+startpos,(kk-startpos)*sizeof(integer));
       
       // upper triangular part used for ETREE0
       // shorten and re-allocate column
       BU.ncol[jj]=startpos;
       BU.rowind[jj]=(integer *)realloc(ja,(size_t)startpos*sizeof(integer));
   } // end for j
#ifdef PRINT_INFO1
   printf("5. pass completed\n");
#endif
   // PRINTMATRIX(&BU);
   // PRINTMATRIX(&BL);

   n=BU.nr;
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
   printf("computed compressed supernode sizes\n");
   for (l=0; l<*nblocks; l++)
       printf("%6ld",blocksize[l]);
   printf("\n"); fflush(stdout);
#endif
   
   SPARSEDELETE(&BL);

   // now expand block information to the original system
   // switch from block size to counter
   idx[0]=0;
   for (l=0; l<*nblocks; l++) 
       idx[l+1]=idx[l]+blocksize[l];
   for (l=0; l<*nblocks; l++) {
       // start position block l, this refers to an index ii in {0,...,BL.nr-1}
       ii=idx[l];
       // this refers to the original index i in {0,...,A->nr-1}
       i=blockstructure[ii];
       // start position block l+1, this refers to an index jj in {1,...,BL.nr}
       jj=idx[l+1];
       // this refers to the original index j in {1,...,A->nr}
       j=blockstructure[jj];
       // their distance is the block size in the original system
       blocksize[l]=j-i;
   } // end for l
#ifdef PRINT_INFO
   printf("computed original supernode sizes\n");
   for (l=0; l<*nblocks; l++)
       printf("%6ld",blocksize[l]);
   printf("\n"); fflush(stdout);
#endif
   
   
   free(idx);
   free(idxpos);
   free(idxpos2);
   free(parent);
   for (l=0; l<n; l++)
       if (nreachable[l]!=0)
	  free(reachable[l]);
   free(reachable);
   free(nreachable);
   free(blockstructure);
   
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


