/* $Id: symselbinv.c 7315 2021-05-28 21:00:20Z bolle $ 

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	Mai 11, 2020. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2020 by TU Braunschweig.  All Rights Reserved.

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




#ifdef _COMPLEX_SYMMETRIC_

#ifdef _SINGLE_REAL_
#define MYSYMSELBINV         SSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       ssytrf
#define MYSYMTRS       ssytrs
#define MYSYMTRI       ssytri
#define MYSYMCON       ssycon
#elif defined _DOUBLE_REAL_
#define MYSYMSELBINV         DSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       dsytrf
#define MYSYMTRS       dsytrs
#define MYSYMTRI       dsytri
#define MYSYMCON       dsycon
#elif defined _SINGLE_COMPLEX_
#define MYSYMSELBINV         CSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       csytrf
#define MYSYMTRS       csytrs
#define MYSYMTRI       csytri
#define MYSYMCON       csycon
#else // double complex
#define MYSYMSELBINV         ZSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       zsytrf
#define MYSYMTRS       zsytrs
#define MYSYMTRI       zsytri
#define MYSYMCON       zsycon
#endif //-if-elif-else single-real

#else // complex Hermitian

#ifdef _SINGLE_REAL_
#define MYSYMSELBINV         SSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       ssytrf
#define MYSYMTRS       ssytrs
#define MYSYMTRI       ssytri
#define MYSYMCON       ssycon
#elif defined _DOUBLE_REAL_
#define MYSYMSELBINV         DSYMselbinv
#define CONJG(A)       (A)
#define MYSYMTRF       dsytrf
#define MYSYMTRS       dsytrs
#define MYSYMTRI       dsytri
#define MYSYMCON       dsycon
#elif defined _SINGLE_COMPLEX_
#define MYSYMSELBINV         CHERselbinv
#define CONJG(A)       (-(A))
#define MYSYMTRF       chetrf
#define MYSYMTRS       chetrs
#define MYSYMTRI       chetri
#define MYSYMCON       checon
#else // double complex
#define MYSYMSELBINV         ZHERselbinv
#define CONJG(A)       (-(A))
#define MYSYMTRF       zhetrf
#define MYSYMTRS       zhetrs
#define MYSYMTRI       zhetri
#define MYSYMCON       zhecon
#endif //-if-elif-else single-real

#endif //-else complex-symmetric



#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

// relaxation parameters when re-allocating buffers
#define ELBOW_BUFF      1.2
#define BUFF_EXT        100

#define MAX_FIELDS 100
#define ELBOW    MAX(4.0,2.0)

// #define PRINT_CHECK 
// #define PRINT_INFO 
// #define printf mexPrintf
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif


void MYSYMSELBINV(SPARSEMATRIX *iA,
		  SPARSEBLOCKMATRIX *BL, SPARSEBLOCKMATRIX *BiD,
		  SPARSEBLOCKMATRIX *BUT,
		  REALS *SL, REALS *SR, integer *perm, integer *invperm)  
{
  integer              i,j,k,l,m,n=BL->nc,p,q, r,s,t, *ncol,nz,
                       ii,jj,flag, size_gemm_buff, sumn, 
                       level3_BLAS, copy_cnt, *block, nblocks=BL->nblocks,
                       n_size, m_size, ni_size, mi_size,
                       i_first, j_first, k_first,
                       Ji_cont, Ik_cont, Ii_cont;
  integer              *piBLJ,     *piBLI,     *piiAI,
                       *piBLinvJ,  *piBLinvI,  
                       *piBLinvJi, *piBLinvIi;
  FLOAT                val, alpha, beta, *gemm_buff, *pviAv,
                       *pv, *pv2, *pv3, *pv4,
                       *pvBLL,     *pvBLD,
                       *pvBLinvL,  *pvBLinvD,
                       *pvBLinvLi, *pvBLinvDi;
    REALS              deltai, deltaj;
    SPARSEBLOCKMATRIX  BLinv;    
#ifdef _PROFILING_
   double timeBegin,
          timeBeginLocal,
          time_symselbinv,
          time_selbinv=0.0,
          time_conversion;

   time_symselbinv=omp_get_wtime();
#endif


    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // initially we simply compute a selected inverse BLinv ~ BL^{-*} BiD BL^{-1}
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    
    // set up BLinv structure, BLinv ~ BL^{-*} BiD BL^{-1}
    BLinv.nblocks=nblocks;
    BLinv.nr=BLinv.nc=n;
    BLinv.nblockcol=BL->nblockcol;
    BLinv.nblockrow=BL->nblockrow;
    BLinv.colind   =BL->colind;
    BLinv.rowind   =BL->rowind;
    BLinv.valD     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *));  // numerical values of the square diagonal blocks
    BLinv.valE     =(FLOAT **)  malloc((size_t)nblocks*sizeof(FLOAT *)); // numerical values of the sub-diagonal blocks

    // auxiliary buffer for level 3 BLAS
    size_gemm_buff=n;
    gemm_buff=(FLOAT *)  malloc((size_t)size_gemm_buff*sizeof(FLOAT));

    // auxiliary buffer for mapping column index i -> block column j
    block    =(integer *)calloc((size_t)n+1,sizeof(integer));
    for (j=0; j<nblocks; j++) {
	// set up block reference for block column j
	piBLJ=BL->colind[j];
	// diagonal block size of BL{j}
	n_size=BL->nblockcol[j];
	for (k=0; k<n_size; k++) {
	    i=piBLJ[k];
	    // block j covers index i
	    block[i]=j;
#ifdef PRINT_CHECK
	    if (i<0 || i>=n)
	       printf(" 1. scalar index out of range when referring to block[%ld]=%ld\n",i+1,j+1); fflush(stdout);
	    if (j<0 || j>=nblocks)
	       printf(" 1. block index out of range when referring to block[%ld]=%ld\n",i+1,j+1); fflush(stdout);
#endif
	} // end for k
    } // end for j 

#ifdef PRINT_INFO
    printf("MYSYMSELBINV: inverse mapping index -> block number computed\n");fflush(stdout);
    for (i=0; i<n; i++)
        printf("%4d", block[i]+1);
    printf("\n");fflush(stdout);
#endif



    // start selected block inversion from the back
    k=nblocks-1;

    // column indices BLinv{k}/BL{k}
    n_size  =BLinv.nblockcol[k];
    piBLinvJ=BLinv.colind[k];
    // create data for diagonal block
    pvBLinvD=BLinv.valD[k]=(FLOAT *)malloc((size_t)n_size*n_size*sizeof(FLOAT));
    // trailing block does not have a sub-diagonal block
    m_size  =0;
    piBLinvI=NULL;
    pvBLinvL=BLinv.valE[k]=NULL;

    // now copy inverse diagonal block information 
    memcpy(pvBLinvD, BiD->valD[k], (size_t)n_size*n_size*sizeof(FLOAT));

    // successively downdate "n" by the size "n_size" of the diagonal block
    sumn=n-n_size;
    
#ifdef PRINT_INFO
    printf("MYSYMSELBINV: final inverse diagonal block computed\n");fflush(stdout);
    printf("        ");
    for (j=0; j<n_size; j++)
        printf("%8ld",piBLinvJ[j]+1);
    printf("\n");fflush(stdout);
    for (i=0; i<n_size; i++) {
        printf("%8ld",piBLinvJ[i]+1);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	for (j=0; j<n_size; j++)
	    printf("%8.1le", pvBLinvD[i+j*n_size]);
#else
	for (j=0; j<n_size; j++)
	    printf("%8.1le", pvBLinvD[i+j*n_size].r);
	printf("\n");
	for (j=0; j<n_size; j++)
	    printf("%8.1le", pvBLinvD[i+j*n_size].i);
#endif
	printf("\n");
	fflush(stdout);
    } // end for i
#endif


    // advance backwards toward to the top
    k--;

#ifdef _PROFILING_
    timeBegin = omp_get_wtime();
#endif   
    // main loop
    while (k>=0) {

          // extract BL{k}
	  piBLJ=BL->colind[k];
	  piBLI=BL->rowind[k];
	  pvBLL=BL->valE[k];

          // diagonal block
          // column indices BLinv{k}/BL{k}
          n_size  =BLinv.nblockcol[k];
	  piBLinvJ=BLinv.colind[k];
	  // create data for diagonal block
	  pvBLinvD=BLinv.valD[k]=(FLOAT *)malloc((size_t)n_size*n_size*sizeof(FLOAT));

	  // sub-diagonal block
          // row indices BLinv{k}/BL{k}
	  m_size  =BLinv.nblockrow[k];
	  piBLinvI=BLinv.rowind[k];
	  // create empty BLinv{k}.L
	  pvBLinvL=BLinv.valE[k]=(FLOAT *)malloc((size_t)m_size*n_size*sizeof(FLOAT));
	  // init with zeros
	  i=m_size*n_size;
	  for (j=0; j<i; j++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      pvBLinvL[j]=0.0;
#else
	      pvBLinvL[j].r=0.0; pvBLinvL[j].i=0.0;
#endif
	  } // end for j
	  // scan the indices of BL{k}.I to find out which block columns are required
	  l=0;
	  while (l<m_size) {
	        // associated row index 
	        ii=piBLI[l];
		// to which block column i does row index ii belong?
	        i=block[ii];
		
		// find out how many row indices of BL{k}.I are associated with block column i
		j=l+1;
		flag=-1;
		while (flag) {
		      if (j>=m_size) {
			 j=m_size-1;
			 flag=0;
		      }
		      else {
			 // associated row index
			 ii=piBLI[j];
			 // end of block column ii exceeded?
			 if (block[ii]>i) { // yes
			    j--;
			    flag=0;
			 }
			 else // no
			    j++;
		      } // end if-else j>=m_size
		} // end while flag 
		// now the row indices within the range l:j of BL{k} are associated with block column BLinv{i}

		// extract already computed BLinv{i}, i>k
		// column index list BLinv{i}.J
		ni_size  =BLinv.nblockcol[i];
		piBLinvJi=BLinv.colind[i];
		// sub-diagonal row index list BLinv{i}.I
		mi_size  =BLinv.nblockrow[i];
		piBLinvIi=BLinv.rowind[i];
		// BLinv{i}.L
		pvBLinvLi=BLinv.valE[i];
		// BLinv{i}.D
		pvBLinvDi=BLinv.valD[i];

		
		// l:j refers to continuously chosen indices !!!
		// Ji, Ik and Ii may exclude some entries !!!


		// check if I(l:j)==I(l):I(j) (contiguous sequence of indices)
		// flag for contiguous index set
		Ji_cont=-1;
		// BLinv{i}.D(Ji,Ji) will physically start at position j_first, 
		// where Ji refers to the sequence of positions in BLinv{i}.D 
		// associated with I(l:j) 
		j_first=piBLI[l]-piBLinvJi[0];
		for (jj=l; jj<=j; jj++) {
		    // index I[jj]
		    ii=piBLI[jj];
		    // non-contiguous index found, break!
		    if (ii>piBLI[l]+jj-l) {
		       Ji_cont=0;
		       jj=j+1;
		    } // end if
		} // end for jj
#ifdef PRINT_INFO
		if (Ji_cont)
		   printf("BL{%d}.I(%d:%d) is a contiguous subsequence of BLinv{%d}.J\n",
			  k+1,l+1,j+1,i+1);
		else
		   printf("BL{%d}.I(%d:%d) does not refer to a contiguous subsequence of BLinv{%d}.J\n",
			  k+1,l+1,j+1,i+1);
		fflush(stdout);
#endif

		// check if the intersection of BLinv{k}.I and BLinv{i}.I
		// consists of contiguous indices
		Ik_cont=-1; Ii_cont=-1;
		p=0; q=0;
		t=0;
		k_first=0; i_first=0;
		while (p<m_size && q<mi_size) {
		      // indices in C-style
		      ii=piBLI[p];
		      jj=piBLinvIi[q];
		      if (ii<jj) {
			 p++;
			 // If we already have common indices, BLinv{k}.I[p]<BLinv{i}.I[q] refers
			 // to a gap in the intersection w.r.t. BLinv{k}.I
		      }
		      else if (ii>jj) {
			 q++;
		      }
		      else { // indices match
			 // store number of the first common index 
			 if (Ik_cont==-1) {
			    // BLinv{k}.L(Ik,:) will physically start at position
			    // k_first, where Ik refers to the sequence of positions
			    // in BLinv{k}.L associated with the intersection of 
			    // BLinv{k}.I and BLinv{i}.I
			    k_first=p;
			    // BLinv{i}.L(Ii,:) will physically start at position
			    // i_first, where Ii refers to the sequence of positions
			    // in BLinv{i}.L associated with the intersection of 
			    // BLinv{k}.I and BLinv{i}.I
			    i_first=q;
			    // store positions of the next indices to stay contiguous
			    Ik_cont=p+1;
			    Ii_cont=q+1;
			 }
			 else {
			    // there exists at least one common index
			    // check if the current index position is the
			    // successor of the previous position 
			    if (p==Ik_cont)
			       // store position of the next index to stay contiguous
			       Ik_cont=p+1;
			    else 
			       Ik_cont=0;
			    if (q==Ii_cont)
			       // store position of the next index to stay contiguous
			       Ii_cont=q+1;
			    else 
			       Ii_cont=0;
			 }
			 p++; q++; t++;
		      } // end if-elseif-else
		} // end while p&q
#ifdef PRINT_INFO
		printf("BL{%d}.I\n",k+1);
		for (p=0; p<m_size; p++)
		    printf("%4d",piBLI[p]+1);
		printf("\n"); 
		fflush(stdout);
		printf("BLinv{%d}.I\n",i+1);
		for (q=0; q<mi_size; q++)
		    printf("%4d",piBLinvIi[q]+1);
		printf("\n"); 
		fflush(stdout);
		if (Ik_cont)
		   printf("intersection leads to a contiguous sequence inside BL{%d}.I of length %d\n",
			  k+1,t);
		else
		   printf("intersection does not yield a contiguous sequence of BL{%d}.I\n",
			  k+1);
		if (Ii_cont)
		   printf("intersection leads to a contiguous sequence inside BLinv{%d}.I  of length %d\n",
			  i+1,t);
		else
		   printf("intersection does not yield a contiguous sequence of BLinv{%d}.I\n",
			  i+1);
		fflush(stdout);
#endif

	  
		// optimal case, all index sets refer to successively stored rows and columns.
		// We can easily use Level 3 BLAS
		if (Ii_cont && Ik_cont && Ji_cont) {
#ifdef PRINT_INFO
		   printf("ideal case, use level 3 BLAS directly!\n");
		   fflush(stdout);
#endif
		   // contribution from the strict lower triangular part
		   // BLinv{k}.L(Ik,:)  = - BLinv{i}.L(Ii,Ji) *BL{k}.L(l:j,:)  + BLinv{k}.L(Ik,:)
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		   alpha=-1.0; beta=1.0;
#else
		   alpha.r=-1.0; beta.r=1.0;
		   alpha.i= 0.0; beta.i=0.0;
#endif
		   ii=j-l+1;
		   if (t && n_size && ii)
		      GEMM("N","N",&t,&n_size,&ii,
			   &alpha,
			   pvBLinvLi+i_first+mi_size*j_first,&mi_size,
			   pvBLL+l,&m_size,
			   &beta,
			   pvBLinvL+k_first,&m_size,1,1);
#ifdef PRINT_INFO
		   printf("Ik=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", r+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   printf("Ii=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", s+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   printf("Ji=[");
		   r=l; s=0;
		   while (r<m_size && s<ni_size) {
		         if (piBLinvJi[s]==piBLI[r]) { 
			    printf("%8d", s+1);
			    r++;
			 }
			 s++;
		   }
		   printf("];\n");
		   printf("MYSYMSELBINV: BLinv{%d}.L(Ik,:) = - BLinv{%d}.L(Ii,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(Ik,:)\n",k+1,i+1,k+1,l+1,j+1,k+1);
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			    for (jj=0; jj<n_size; jj++)
			        printf("%8.1le",pvBLinvL[r+m_size*jj]);
#else
			    for (jj=0; jj<n_size; jj++)
			        printf("%8.1le",pvBLinvL[r+m_size*jj].r);
			    printf("\n");
			    for (jj=0; jj<n_size; jj++)
			        printf("%8.1le",pvBLinvL[r+m_size*jj].i);
#endif
			    printf("\n");
			    fflush(stdout);
			    r++; s++; 
			 } // end if-elseif-else
		   } // end while r&s
#endif

		   // contribution from the strict upper triangular part 
		   // BLinv{k}.L(l:j,:) = - BLinv{i}.L(Ii,Ji)^T*BL{k}.L(Ik,:)  + BLinv{k}.L(l:j,:)
		   if (ii && n_size && t) {
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_ &&  !defined _COMPLEX_SYMMETRIC_
		      // complex-Hermitian case
		      GEMM("C","N",&ii,&n_size,&t,
			   &alpha,
			   pvBLinvLi+i_first+mi_size*j_first,&mi_size,
			   pvBLL+k_first,&m_size,
			   &beta,
			   pvBLinvL+l,&m_size,1,1);
#else 
		      GEMM("T","N",&ii,&n_size,&t,
			   &alpha,
			   pvBLinvLi+i_first+mi_size*j_first,&mi_size,
			   pvBLL+k_first,&m_size,
			   &beta,
			   pvBLinvL+l,&m_size,1,1);
#endif
		   } // end if
#ifdef PRINT_INFO
		   printf("Ik=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", r+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   printf("Ii=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", s+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   printf("Ji=[");
		   r=l; s=0;
		   while (r<m_size && s<ni_size) {
		         if (piBLinvJi[s]==piBLI[r]) { 
			    printf("%8d", s+1);
			    r++;
			 }
			 s++;
		   }
		   printf("];\n");
		   printf("MYSYMSELBINV: BLinv{%d}.L(%d:%d,:) = - BLinv{%d}.L(Ii,Ji)' *BL{%d}.L(Ik,:)  + BLinv{%d}.L(%d:%d,:)\n",k+1,l+1,j+1,i+1,k+1,k+1,l+1,j+1);
		   for (jj=l; jj<=j; jj++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		       for (q=0; q<n_size; q++)
			   printf("%8.1le",pvBLinvL[jj+m_size*q]);
#else
		       for (q=0; q<n_size; q++)
			   printf("%8.1le",pvBLinvL[jj+m_size*q].r);
		       printf("\n");
		       for (q=0; q<n_size; q++)
			   printf("%8.1le",pvBLinvL[jj+m_size*q].i);
#endif
		       printf("\n");
		       fflush(stdout);
		   }
#endif

		   // contribution from the diagonal block
		   // BLinv{k}.L(l:j,:) = - BLinv{i}.D(Ji,Ji)  *BL{k}.L(l:j,:) + BLinv{k}.L(l:j,:)
		   if (ii && n_size)
		      GEMM("N","N",&ii,&n_size,&ii,
			   &alpha,
			   pvBLinvDi+j_first+j_first*ni_size,&ni_size,
			   pvBLL+l,&m_size,
			   &beta,
			   pvBLinvL+l,&m_size,1,1);
#ifdef PRINT_INFO
		   r=l; s=0;
		   printf("Ji=[");
		   while (r<m_size && s<ni_size) {
		         if (piBLinvJi[s]==piBLI[r]) { 
			    printf("%8d", s+1);
			    r++;
			 }
			 s++;
		   }
		   printf("];\n");
		   printf("MYSYMSELBINV: BLinv{%d}.L(%d:%d,:) = - BLinv{%d}.D(Ji,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(%d:%d,:)\n",
			     k+1,l+1,j+1,i+1,k+1,l+1,j+1,k+1,l+1,j+1);
		   for (r=l; r<=j; r++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s]);
#else
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s].r);
		       printf("\n");
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s].i);
#endif
		       printf("\n");
		       fflush(stdout);
		   }
#endif
		} // end if Ii_cont & Ik_cont & Ji_cont
		else { // now at least one block is not contiguous. The decision 
		       // whether to stik with level 3 BLAS or not will be made on
		       // the cost for copying part of the data versus the 
		       // computational cost. This is definitely not optimal
		   


		   // ----------------------------------------------------------
  		   // ----------------------------------------------------------
		   // --- contribution from the strict lower triangular part ---
		   // BLinv{k}.L(Ik,:)  = - BLinv{i}.L(Ii,Ji)  *BL{k}.L(l:j,:) + BLinv{k}.L(Ik,:)
		   // determine amount of auxiliary memory
#ifdef PRINT_INFO
		   if (!Ji_cont)
		      printf("Ji not contiguous\n");
		   if (!Ii_cont)
		      printf("Ii not contiguous\n");
		   if (!Ik_cont)
		      printf("Ik not contiguous\n");
		   fflush(stdout);
#endif
		   copy_cnt=0;
		   // level 3 BLAS have to use |Ii| x |Ji| buffer rather than BLinv{i}.L(Ii,Ji)
		   if (!Ii_cont || !Ji_cont)
		      copy_cnt+=t*(j-l+1);
		   // level 3 BLAS have to use |Ik| x n_size buffer rather than BLinv{k}.L(Ik,:)
		   if (!Ik_cont)
		      copy_cnt+=t*n_size;

		   if (copy_cnt<t*(j-l+1)*n_size)
		      level3_BLAS=-1;
		   else 
		      level3_BLAS=0;
		
		   // it could pay off to copy the data into one or two auxiliary buffers
		   if (level3_BLAS && t) {
#ifdef PRINT_INFO
		      printf("contribution from the strict lower triangular part, still use level 3 BLAS\n");
		      fflush(stdout);
#endif
		      size_gemm_buff=MAX(size_gemm_buff,copy_cnt);
		      gemm_buff=(FLOAT *)realloc(gemm_buff,
						 (size_t)size_gemm_buff*sizeof(FLOAT));
		      if (!Ii_cont || !Ji_cont) {
			 // copy BLinv{i}.L(Ii,Ji) to buffer
		         pv=gemm_buff;
			 p=0; q=0;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
			 	  p++;
			       else if (ii>jj)
				  q++;
			       else { // indices match

				  // copy parts of the current row BLinv{i}.L(p,:) of 
				  // BLinv{i}.L(Ii,Ji) associated with Ji to gemm_buff
				  pv3=pv;
				  pv2=pvBLinvLi+q;
				  r=l; s=0;
				  while (r<m_size && s<ni_size) {
				        // does column BL{k}.I(r) match some BLinv{i}.J(s)?
					// Recall that I(l:j) is a subset of Ji 
				        if (piBLinvJi[s]==piBLI[r]) { 
					   *pv3=*pv2;
					   pv3+=t;
					   r++;
					} // end if
					s++;
					pv2+=mi_size;
				  } // end while s
				  pv++;
				
				  p++; q++; 
			       } // end if-elseif-else
			 } // end while p&q
#ifdef PRINT_INFO
			 printf("Ik=[");
			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
		               if (piBLI[r]<piBLinvIi[s]) 
			          r++;
			       else if (piBLI[r]>piBLinvIi[s])
				  s++;
			       else { 
		                  printf("%8d", r+1);
				  r++; s++; 
			       }
			 }
			 printf("];\n");
			 printf("Ii=[");
			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
			       if (piBLI[r]<piBLinvIi[s]) 
			          r++;
			       else if (piBLI[r]>piBLinvIi[s])
				  s++;
			       else { 
		                  printf("%8d", s+1);
				  r++; s++; 
			       }
			 }
			 printf("];\n");
			 printf("Ji=[");
			 r=l; s=0;
			 while (r<m_size && s<ni_size) {
		               if (piBLinvJi[s]==piBLI[r]) { 
		                  printf("%8d", s+1);
				  r++;
			       }
			       s++;
			 }
			 printf("];\n");
			 printf("MYSYMSELBINV: BLinv{%d}.L(Ii,Ji) cached\n",i+1);fflush(stdout);
			 printf("        ");
			 r=l; s=0;
			 while (r<m_size && s<ni_size) {
			       if (piBLinvJi[s]==piBLI[r]) {
				  printf("%8d", piBLinvJi[s]+1);
				  r++;
			       }
			       s++;
			 } // end while s
			 printf("\n");fflush(stdout);
			 p=0; q=0;
			 pv=gemm_buff;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
				  p++;
			       else if (ii>jj)
				  q++;
			       else { // indices match
				  printf("%8d", ii+1);

				  r=l; s=0;
				  pv2=pv;
				  while (r<m_size && s<ni_size) {
				        if (piBLinvJi[s]==piBLI[r]) { 
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
					   printf("%8.1le", *pv2);
#else
					   printf("%8.1le", pv2->r);
					   printf("+%8.1lei", pv2->i);
#endif
					   pv2+=t;
					   r++;
					}
					s++;
				  }
				  pv++;
				  printf("\n");fflush(stdout);

				  p++; q++; 
			       }
			 }
#endif

			 pv=gemm_buff; p=t;
		      } // end if !Ii_cont or !Ji_cont
		      else { 
			 // pointer to BLinv{i}.L(Ii,Ji) and LDA
		         pv=pvBLinvLi+i_first+mi_size*j_first; p=mi_size;
		      } // end if-else
		   
		      if (!Ik_cont) {
		         // init buffer with zeros
		         if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;
			 for (q=0; q<t*n_size; q++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			     *pv2++=0.0;
#else
			     pv2->r=0.0;
			     pv2->i=0.0;
			     pv2++;
#endif
			 } // end for q
			 // pointer and LDC
			 if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;
			 q=t;
			 // since we initialized everything with zero, beta is 
			 // almost arbitrary, we indicate this changing beta to 0 
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			 alpha=1.0; beta=0.0;
#else
			 alpha.r=1.0; beta.r=0.0;
			 alpha.i=0.0; beta.i=0.0;
#endif
#ifdef PRINT_INFO
			 printf("MYSYMSELBINV: cached zeros instead of  BLinv{%d}.L(Ik,:)\n",k+1);
			 fflush(stdout);
#endif
		      } // end if !Ik_cont
		      else {
			 // pointer to BLinv{k}.L(Ik,:) and LDC
		         pv2=pvBLinvL+k_first; q=m_size;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			 alpha=-1.0; beta=1.0;
#else
			 alpha.r=-1.0; beta.r=1.0;
			 alpha.i= 0.0; beta.i=0.0;
#endif
		      } // end if-else !Ik_cont
		      
                      // call level 3 BLAS
		      ii=j-l+1;
		      if (t && n_size && ii)
			 GEMM("N","N",&t,&n_size,&ii,
			      &alpha,
			      pv,&p,
			      pvBLL+l,&m_size,
			      &beta,
			      pv2,&q,1,1);
#ifdef PRINT_INFO
		      printf("Ik=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
		            if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", r+1);
			       r++; s++; 
		            }
		      }
		      printf("];\n");
		      printf("Ii=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
			    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", s+1);
			       r++; s++; 
		            }
		      }
		      printf("];\n");
		      r=l; s=0;
		      printf("Ji=[");
		      while (r<m_size && s<ni_size) {
		            if (piBLinvJi[s]==piBLI[r]) { 
		               printf("%8d", s+1);
			       r++;
		            }
			    s++;
		      }
		      printf("];\n");
		      if (Ik_cont)
			 printf("MYSYMSELBINV: BLinv{%d}.L(Ik,:) = - BLinv{%d}.L(Ii,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(Ik,:)\n",
		                   k+1,i+1,k+1,l+1,j+1,k+1);
		      else
			 printf("MYSYMSELBINV: cached                BLinv{%d}.L(Ii,Ji)  *BL{%d}.L(%d:%d,:)\n",
		                   i+1,k+1,l+1,j+1);

		      for (r=0; r<t; r++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		          for (s=0; s<n_size; s++)
			      printf("%8.1le",pv2[r+q*s]);
#else
		          for (s=0; s<n_size; s++)
			      printf("%8.1le",pv2[r+q*s].r);
			  printf("\n");
		          for (s=0; s<n_size; s++)
			      printf("%8.1le",pv2[r+q*s].i);
#endif
			  printf("\n");
			  fflush(stdout);
		      }
#endif
		     
		     
		      if (!Ik_cont) {
			 // init buffer with zeros
		         if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;
			 p=0; q=0;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
				  p++;
			       else if (ii>jj)
				  q++;
			       else { // indices match

				  // copy current row of pv2 to BLinv{k}.L(Ik,:)
				  pv=pvBLinvL+p;
				  pv3=pv2;
				  for (r=0; r<n_size; r++, pv+=m_size, pv3+=t) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
				      *pv-=*pv3;
#else
				      pv->r-=pv3->r;
				      pv->i-=pv3->i;
#endif
				  } // end for r
				  pv2++;
				
				  p++; q++; 
			       } // end if-elseif-else
			 } // end while p&q
#ifdef PVINT_INFO
			 printf("Ik=[");
			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
			       if (piBLI[r]<piBLinvIi[s]) 
			          r++;
			       else if (piBLI[r]>piBLinvIi[s])
			          s++;
			       else { 
		                  printf("%8d", r+1);
			          r++; s++; 
		               }
		         }
			 printf("];\n");
			 printf("Ii=[");
			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
			       if (piBLI[r]<piBLinvIi[s]) 
			          r++;
			       else if (piBLI[r]>piBLinvIi[s])
			          s++;
			       else { 
		                  printf("%8d", s+1);
				  r++; s++; 
		               }
		         }
			 printf("];\n");
			 printf("Ji=[");
			 r=l; s=0;
			 while (r<m_size && s<ni_size) {
		               if (piBLinvJi[s]==piBLI[r]) { 
		                  printf("%8d", s+1);
				  r++;
		               }
			       s++;
		         }
			 printf("];\n");
			 printf("MYSYMSELBINV: BLinv{%d}.L(Ik,:) = - BLinv{%d}.L(Ii,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(Ik,:)\n",k+1,i+1,k+1,l+1,j+1,k+1);

			 p=0; q=0;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
			 	  p++;
			       else if (ii>jj)
				  q++;
			       else {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		                  for (s=0; s<n_size; s++)
			              printf("%8.1le",pvBLinvL[p+m_size*s]);
#else
		                  for (s=0; s<n_size; s++)
			              printf("%8.1le",pvBLinvL[p+m_size*s].r);
				  printf("\n");
		                  for (s=0; s<n_size; s++)
			              printf("%8.1le",pvBLinvL[p+m_size*s].i);
#endif
				  printf("\n");
				  fflush(stdout);
				  p++; q++; 
		               } // end if-elseif-else
		         } // end while p&q
#endif


		      } // if !Ik_cont
		   } // end if level3_BLAS
		   else if (t) { // it might not pay off, therefore we use a simple hand-coded loop
		      // BLinv{k}.L(Ik,:)  -=  BLinv{i}.L(Ii,Ji) * BL{k}.L(l:j,:)
#ifdef PRINT_INFO
		      printf("contribution from the strict lower triangular part, use hand-coded loops\n");
		      fflush(stdout);
#endif
		      p=0; q=0;
		      while (p<m_size && q<mi_size) {
			    ii=piBLI[p];
			    jj=piBLinvIi[q];
			    if (ii<jj)
			       p++;
			    else if (ii>jj)
			       q++;
			    else { // row indices BLinv{k}.I[p]=BLinv{i}.I[q] match
			       pv =pvBLinvL+p;
			       pv3=pvBLL+l;
			       // BLinv{k}.L(p,:)  -=  BLinv{i}.L(q,Ji) * BL{k}.L(l:j,:)
			       for (ii=0; ii<n_size; ii++,pv+=m_size,pv3+=m_size-(j-l+1)) {
				   // BLinv{k}.L(p,ii)  -=  BLinv{i}.L(q,Ji) * BL{k}.L(l:j,ii)
				   pv2=pvBLinvLi+q;
				   r=l; s=0;
				   while (r<m_size && s<ni_size) {
				         // column Ji[s] occurs within I(l:j). 
					 // Recall that I(l:j) is a subset of Ji 
				         if (piBLinvJi[s]==piBLI[r]) { 
					    // BLinv{k}.L(p,ii)  -=  BLinv{i}.L(q,s) * BL{k}.L(r,ii)
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
					    *pv-=(*pv2)*(*pv3++);
#else
					    val.r=pv2->r*pv3->r-pv2->i*pv3->i;
					    val.i=pv2->r*pv3->i+pv2->i*pv3->r;
					    pv->r-=val.r;
					    pv->i-=val.i;
					    pv3++;
#endif
					    r++;
					 }
					 s++;
					 pv2+=mi_size;
				   } // end while s
			       } // end for ii
			       p++; q++; 
			    } // end if-elseif-else
		      } // end while p&q
#ifdef PRINT_INFO
		      printf("Ik=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
			    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
			       printf("%8d", r+1);
			       r++; s++; 
			    }
		      }
		      printf("];\n");
		      printf("Ii=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
			    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
			       printf("%8d", s+1);
			       r++; s++; 
			    }
		      }
		      printf("];\n");
		      printf("Ji=[");
		      r=l; s=0;
		      while (r<m_size && s<ni_size) {
			    if (piBLinvJi[s]==piBLI[r]) { 
			       printf("%8d", s+1);
			       r++;
			    }
			    s++;
		      }
		      printf("];\n");
		      printf("MYSYMSELBINV: BLinv{%d}.L(Ik,:) = - BLinv{%d}.L(Ii,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(Ik,:)\n",k+1,i+1,k+1,l+1,j+1,k+1);
		      p=0; q=0;
		      while (p<m_size && q<mi_size) {
			    ii=piBLI[p];
			    jj=piBLinvIi[q];
			    if (ii<jj) 
			       p++;
			    else if (ii>jj)
			       q++;
			    else {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			       for (s=0; s<n_size; s++)
				   printf("%8.1le",pvBLinvL[p+m_size*s]);
#else
			       for (s=0; s<n_size; s++)
				   printf("%8.1le",pvBLinvL[p+m_size*s].r);
			       printf("\n");
			       for (s=0; s<n_size; s++)
				   printf("%8.1le",pvBLinvL[p+m_size*s].i);
#endif
			       printf("\n");
			       fflush(stdout);
			       p++; q++; 
			    } // end if-elseif-else
		      } // end while p&q
#endif
		   } // end if-else level3_BLAS
		   else {
#ifdef PRINT_INFO
		      printf("contribution from the strict lower triangular part empty\n");
		      fflush(stdout);
#endif
		   }
		   // end contribution from the strict lower triangular part
		   // ------------------------------------------------------
		   // ------------------------------------------------------

		

		   // ----------------------------------------------------------
		   // ----------------------------------------------------------
		   // --- contribution from the strict upper triangular part ---
		   // BLinv{k}.L(l:j,:) = - BLinv{i}.L(Ii,Ji)^T*BL{k}.L(Ik,:)  + BLinv{k}.L(l:j,:)
		   // determine amount of auxiliary memory
		   copy_cnt=0;
		   // level 3 BLAS have to use |Ii| x |Ji| buffer rather than BLinv{i}.L(Ii,Ji)
		   if (!Ii_cont || !Ji_cont)
		      copy_cnt+=t*(j-l+1);
		   // level 3 BLAS have to use |Ik| x n_size buffer rather than BLinv{k}.L(Ik,:)
		   if (!Ik_cont)
		      copy_cnt+=t*n_size;

		   if (copy_cnt<t*(j-l+1)*n_size)
		      level3_BLAS=-1;
		   else 
		      level3_BLAS=0;

		   // it could pay off to copy the data into one or two auxiliary buffers
		   if (level3_BLAS && t) {
#ifdef PRINT_INFO
		      printf("contribution from the strict upper triangular part, still use level 3 BLAS\n");
		      fflush(stdout);
#endif
		      size_gemm_buff=MAX(size_gemm_buff,copy_cnt);
		      gemm_buff=(FLOAT *)realloc(gemm_buff,
						 (size_t)size_gemm_buff*sizeof(FLOAT));
		      if (!Ii_cont || !Ji_cont) {
			 // copy BLinv{i}.L(Ii,Ji) to buffer
		         pv=gemm_buff;
			 p=0; q=0;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
			 	  p++;
			       else if (ii>jj)
				  q++;
			       else { // indices match

				  // copy parts of the current row of BLinv{i}.L(Ii,:) 
				  // associated with Ji to gemm_buff
				  pv3=pv;
				  pv2=pvBLinvLi+q;
				  r=l; s=0;
				  while (r<m_size && s<ni_size) {
				        // column Ji[s] occurs within I(l:j). 
					// Recall that I(l:j) is a subset of Ji 
				        if (piBLinvJi[s]==piBLI[r]) { 
					   *pv3=*pv2;
					   pv3+=t;
					   r++;
					}
					s++;
					pv2+=mi_size;
				  } // end while s
				  pv++;
				
				  p++; q++; 
			       } // end if-elseif-else
			 } // end while p&q
#ifdef PRINT_INFO
			 printf("MYSYMSELBINV: cached copy of BLinv{%d}.L(Ii,Ji)\nIndex set Ji:\n",i+1);
			 fflush(stdout);
			 r=l; s=0;
			 while (r<m_size && s<ni_size) {
			       if (piBLinvJi[s]==piBLI[r]) { 
				  printf("%8d",piBLinvJi[s]+1);
				  r++;
			       }
			       s++;
			 } // end while s
			 printf("\nIndex set Ii:\n");
			 fflush(stdout);
			 p=0; q=0;
			 while (p<m_size && q<mi_size) {
			       ii=piBLI[p];
			       jj=piBLinvIi[q];
			       if (ii<jj) 
			 	  p++;
			       else if (ii>jj)
				  q++;
			       else {
				  printf("%8d",ii+1);
				  p++; q++; 
			       } // end if-elseif-else
			 } // end while p&q
			 printf("\n");
			 fflush(stdout);
			 for (p=0; p<t; p++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			     for (q=0; q<j-l+1; q++)
			         printf("%8.1le",gemm_buff[p+q*t]);
#else
			     for (q=0; q<j-l+1; q++)
			         printf("%8.1le",gemm_buff[p+q*t].r);
			     printf("\n");
			     for (q=0; q<j-l+1; q++)
			         printf("%8.1le",gemm_buff[p+q*t].i);
#endif
			     printf("\n");
			     fflush(stdout);
			 }
#endif

			 pv=gemm_buff; p=t;
		      } 
		      else {
			 // pointer to BLinv{i}.L(Ii,Ji) and LDA
		         pv=pvBLinvLi+i_first+mi_size*j_first; p=mi_size;
		      } // end if-else

		      if (!Ik_cont) {
			 // copy BL{k}.L(Ik,:) to buffer
		         if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;

			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
			       ii=piBLI[r];
			       jj=piBLinvIi[s];
			       if (ii<jj) 
				  r++;
			       else if (ii>jj)
				  s++;
			       else { // indices match
				 
				  // copy BL{k}.L(r,:) to buffer
				  pv3=pv2;
				  pv4=pvBLL+r;
				  for (ii=0; ii<n_size; ii++,pv3+=t,pv4+=m_size) 
				      *pv3=*pv4;
				  pv2++;
				
				  r++; s++; 
			       } // end if-elseif-else
			 } // end while p&q
#ifdef PRINT_INFO
			 printf("MYSYMSELBINV: cached copy of BL{%d}.L(Ik,:)\nIndex set J:\n",i+1);
			 fflush(stdout);
			 for (q=0; q<n_size; q++)
			     printf("%8d",piBLJ[q]+1);
			 printf("\nIndex set Ik:\n");
			 fflush(stdout);
			 r=0; s=0;
			 while (r<m_size && s<mi_size) {
			       ii=piBLI[r];
			       jj=piBLinvIi[s];
			       if (ii<jj) 
			 	  r++;
			       else if (ii>jj)
				  s++;
			       else {
				  printf("%8d",ii+1);
				  r++; s++; 
			       } // end if-elseif-else
			 } // end while p&q
			 printf("\n");
			 fflush(stdout);
			 if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;
			 for (r=0; r<t; r++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			     for (s=0; s<n_size; s++)
			         printf("%8.1le",pv2[r+s*t]);
#else
			     for (s=0; s<n_size; s++)
			         printf("%8.1le",pv2[r+s*t].r);
			     printf("\n");
			     for (s=0; s<n_size; s++)
			         printf("%8.1le",pv2[r+s*t].i);
#endif
			     printf("\n");
			     fflush(stdout);
			 }
#endif


			 // pointer and LDC
			 if (!Ii_cont || !Ji_cont)
			    pv2=gemm_buff+t*(j-l+1);
			 else
			    pv2=gemm_buff;
			 q=t;
		      } // end if !Ik_cont
		      else {
			 // pointer to BL{k}.L(Ik,:) and LDC
		         pv2=pvBLL+k_first; q=m_size;
		      } // end if-else !Ik_cont

		      // call level 3 BLAS
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		      alpha=-1.0; beta=1.0;
#else
		      alpha.r=-1.0; beta.r=1.0;
		      alpha.i= 0.0; beta.i=0.0;
#endif
		      ii=j-l+1;
		      if (ii && n_size && t)
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_ &&  !defined _COMPLEX_SYMMETRIC_
			 // complex-Hermitian case
			 GEMM("C","N",&ii,&n_size,&t,
			      &alpha,
			      pv,&p,
			      pv2,&q,
			      &beta,
			      pvBLinvL+l,&m_size,1,1);
#else
			 GEMM("T","N",&ii,&n_size,&t,
			      &alpha,
			      pv,&p,
			      pv2,&q,
			      &beta,
			      pvBLinvL+l,&m_size,1,1);
#endif
#ifdef PRINT_INFO
		      printf("Ik=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
		 	    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", r+1);
			       r++; s++; 
		            }
		      }
		      printf("];\n");
		      printf("Ii=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
		 	    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", s+1);
			       r++; s++; 
		            }
		      }
		      printf("];\n");
		      r=l; s=0;
		      printf("Ji=[");
		      while (r<m_size && s<ni_size) {
		            if (piBLinvJi[s]==piBLI[r]) { 
		               printf("%8d", s+1);
			       r++;
		            }
			    s++;
		      }
		      printf("];\n");
		      printf("MYSYMSELBINV: BLinv{%d}.L(%d:%d,:) = - BLinv{%d}.L(Ii,Ji)^T*BL{%d}.L(Ik,:)  + BLinv{%d}.L(%d:%d,:)\n",k+1,l+1,j+1,i+1,k+1,k+1,l+1,j+1);
		      for (jj=l; jj<=j; jj++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[jj+m_size*q]);
#else
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[jj+m_size*q].r);
			  printf("\n");
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[jj+m_size*q].i);
#endif
			  printf("\n");
			  fflush(stdout);
		      }
#endif
		      
		   } // end if level3_BLAS 
		   else if (t) { /* it might not pay off, therefore we use a simple hand-coded loop */
		      // BLinv{k}.L(l:j,:) -=  BLinv{i}.L(Ii,Ji)^T * BL{k}.L(Ik,:)
#ifdef PRINT_INFO
		      printf("contribution from the strict upper triangular part, use hand-coded loops\n");
		      fflush(stdout);
#endif
		      p=0; q=0;
		      while (p<m_size && q<mi_size) {
			    ii=piBLI[p];
			    jj=piBLinvIi[q];
			    if (ii<jj)
			       p++;
			    else if (ii>jj)
			       q++;
			    else { // row indices BL{k}.I[p]=BLinv{i}.I[q] match
			       pv =pvBLL+p;
			       pv3=pvBLinvL+l;
			       // BLinv{k}.L(l:j,:) -=  BLinv{i}.L(q,Ji)^T * BL{k}.L(p,:)
			       for (ii=0; ii<n_size; ii++,pv+=m_size,pv3+=m_size-(j-l+1)) {
				   // BLinv{k}.L(l:j,ii) -=  BLinv{i}.L(q,Ji)^T * BL{k}.L(p,ii)
				   pv2=pvBLinvLi+q;
				   r=l; s=0;
				   while (r<m_size && s<ni_size) {
				         // column Ji[s] occurs within I(l:j). 
					 // Recall that I(l:j) is a subset of Ji 
				         if (piBLinvJi[s]==piBLI[r]) { 
					    // BLinv{k}.L(r,ii)  -=  BLinv{i}.L(q,s)^T * BL{k}.L(p,ii)
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
					    *pv3++ -= (*pv2)*(*pv);
#else
					    val.r=pv2->r*pv->r-pv2->i*pv->i;
					    val.i=pv2->r*pv->i+pv2->i*pv->r;
					    pv3->r-=val.r;
					    pv3->i-=val.i;
					    pv3++;
#endif
					    r++;
					 }
					 s++;
					 pv2+=mi_size;
				   } // end while s
			       } // end for ii
			       p++; q++; 
			    } // end if-elseif-else
		      } // end while p&q
#ifdef PRINT_INFO
		      printf("Ik=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
			    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", r+1);
			       r++; s++; 
			    }
		      }
		      printf("];\n");
		      printf("Ii=[");
		      r=0; s=0;
		      while (r<m_size && s<mi_size) {
			    if (piBLI[r]<piBLinvIi[s]) 
			       r++;
			    else if (piBLI[r]>piBLinvIi[s])
			       s++;
			    else { 
		               printf("%8d", s+1);
			       r++; s++; 
			    }
		      }
		      printf("];\n");
		      r=l; s=0;
		      printf("Ji=[");
		      while (r<m_size && s<ni_size) {
		            if (piBLinvJi[s]==piBLI[r]) { 
		               printf("%8d", s+1);
			       r++;
			    }
			    s++;
		      }
		      printf("];\n");
		      printf("MYSYMSELBINV: BLinv{%d}.L(%d:%d,:) = - BLinv{%d}.L(Ii,Ji)^T*BL{%d}.L(Ik,:)  + BLinv{%d}.L(%d:%d,:)\n",k+1,l+1,j+1,i+1,k+1,k+1,l+1,j+1);
		      for (p=l; p<=j; p++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[p+m_size*q]);
#else
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[p+m_size*q].r);
			  printf("\n");
			  for (q=0; q<n_size; q++)
			      printf("%8.1le",pvBLinvL[p+m_size*q].i);
#endif
			  printf("\n");
			  fflush(stdout);
		      }
#endif
		   } // end if-else level3_BLAS
		   else {
#ifdef PRINT_INFO
		      printf("contribution from the strict upper triangular part empty\n");
		      fflush(stdout);
#endif
		   }
		   // end contribution from the strict upper triangular part
		   // --------------------------------------------------------
		   // --------------------------------------------------------



		   // --------------------------------------------------------
		   // --------------------------------------------------------
		   // ---------  contribution from the diagonal block --------
		   // BLinv{k}.L(l:j,:) = - BLinv{i}.D(Ji,Ji)  *BL{k}.L(l:j,:) + BLinv{k}.L(l:j,:)
		   // determine amount of auxiliary memory
		   copy_cnt=0;
		   // level 3 BLAS have to use |Ji| x |Ji| buffer rather than BLinv{i}.D(Ji,Ji)
		   if (!Ji_cont)
		      copy_cnt=(j-l+1)*(j-l+1);

		   // it pays off to copy the data into one or two auxiliary buffers
		   size_gemm_buff=MAX(size_gemm_buff,copy_cnt);
		   gemm_buff=(FLOAT *)realloc(gemm_buff,(size_t)size_gemm_buff*sizeof(FLOAT));
		   if (!Ji_cont) {
		      // copy BLinv{i}.D(Ji,Ji) to buffer
		      pv=gemm_buff;
		      pv2=pvBLinvDi;
		      p=l; q=0;
		      while (p<m_size && q<ni_size) {
			    // column Ji[q] occurs within I(l:j). 
			    // Recall that I(l:j) is a subset of Ji 
			    if (piBLinvJi[q]==piBLI[p]) { 
			       // copy BLinv{i}.D(Ji,q) to buffer
			       r=l; s=0;
			       while (r<m_size && s<ni_size) {
				     // column Ji[s] occurs within I(l:j). 
				     // Recall that I(l:j) is a subset of Ji 
				     if (piBLinvJi[s]==piBLI[r]) { 
				        *pv++=*pv2;
					r++;
				     }
				     pv2++;

				     s++;
			       } // end while r&s
			       
			       p++; 
			    } // end if
			    else
			       pv2+=ni_size;
			    q++; 
		      } // end while p&q
#ifdef PRINT_INFO
		      printf("MYSYMSELBINV: cached copy of BLinv{%d}.D(Ji,Ji)\nIndices:\n",i+1);
		      fflush(stdout);
		      p=l; q=0;
		      while (p<m_size && q<ni_size) {
			    if (piBLinvJi[q]==piBLI[p]) { 
			       printf("%8d",piBLinvJi[q]+1);
			       p++; 
			    } // end if
			    q++; 
		      } // end while p&q
		      printf("\n");
		      fflush(stdout);
		      for (p=0; p<j-l+1; p++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
			  for (q=0; q<j-l+1; q++)
			      printf("%8.1le",gemm_buff[p+q*(j-l+1)]);
#else
			  for (q=0; q<j-l+1; q++)
			      printf("%8.1le",gemm_buff[p+q*(j-l+1)].r);
			  printf("\n");
			  for (q=0; q<j-l+1; q++)
			      printf("%8.1le",gemm_buff[p+q*(j-l+1)].i);
#endif
			  printf("\n");
			  fflush(stdout);
		      }

#endif
		      pv=gemm_buff; p=j-l+1;
		   } // end if !Ji_cont
		   else {
		      // pointer to BLinv{i}.D(Ji,Ji) and LDA
		      pv=pvBLinvDi+j_first+j_first*ni_size; p=ni_size;
		   } // end if-else  !Ji_cont

		   
		   // call level 3 BLAS
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		   alpha=-1.0; beta=1.0;
#else
		   alpha.r=-1.0; beta.r=1.0;
		   alpha.i= 0.0; beta.i=0.0;
#endif
		   ii=j-l+1;
		   if (ii && n_size)
		      GEMM("N","N",&ii,&n_size,&ii,
			   &alpha,
			   pv,&p,
			   pvBLL+l,&m_size,
			   &beta,
			   pvBLinvL+l,&m_size,1,1);
#ifdef PRINT_INFO
		   printf("Ik=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", r+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   printf("Ii=[");
		   r=0; s=0;
		   while (r<m_size && s<mi_size) {
		         if (piBLI[r]<piBLinvIi[s]) 
			    r++;
			 else if (piBLI[r]>piBLinvIi[s])
			    s++;
			 else { 
			    printf("%8d", s+1);
			    r++; s++; 
			 }
		   }
		   printf("];\n");
		   r=l; s=0;
		   printf("Ji=[");
		   while (r<m_size && s<ni_size) {
		         if (piBLinvJi[s]==piBLI[r]) { 
			    printf("%8d", s+1);
			    r++;
			 }
			 s++;
		   }
		   printf("];\n");
		   printf("MYSYMSELBINV: BLinv{%d}.L(%d:%d,:) = - BLinv{%d}.D(Ji,Ji)  *BL{%d}.L(%d:%d,:)  + BLinv{%d}.L(Ik,:)\n",
			     k+1,l+1,j+1,i+1,k+1,l+1,j+1,k+1);
		   for (r=l; r<=j; r++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s]);
#else
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s].r);
		       printf("\n");
		       for (s=0; s<n_size; s++)
			   printf("%8.1le",pvBLinvL[r+m_size*s].i);
#endif
		       printf("\n");
		       fflush(stdout);
		   }
#endif
		   
		   // end contribution from the strict lower triangular part
		   // ---------------------------------------------------------
		   // ---------------------------------------------------------
		} // end if-else Ii_cont & Ik_cont & Ji_cont
		   
		// advance to the next block column
		l=j+1;
	  } // end while l<m_size

#ifdef PRINT_INFO
	  printf("MYSYMSELBINV: %d-th inverse sub-diagonal block computed\n", k+1);fflush(stdout);
	  printf("        ");
	  for (j=0; j<n_size; j++)
	      printf("%8d", piBLinvJ[j]+1);
	  printf("\n");fflush(stdout);
	  for (i=0; i<m_size; i++) {
	      printf("%8d", piBLinvI[i]+1);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvL[i+j*m_size]);
#else
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvL[i+j*m_size].r);
	      printf("\n");
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvL[i+j*m_size].i);
#endif
	      printf("\n");
	      fflush(stdout);
	  }
#endif

	  // now copy inverse diagonal block information 
	  memcpy(pvBLinvD, BiD->valD[k], (size_t)n_size*n_size*sizeof(FLOAT));



	  // BLinv{k}.D = - BL{k}.L^T *BLinv{k}.L + BLinv{k}.D
	  // call level 3 BLAS
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  alpha=-1.0; beta=1.0;
#else
	  alpha.r=-1.0; beta.r=1.0;
	  alpha.i= 0.0; beta.i=0.0;
#endif
	  if (n_size && m_size) {
#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_ && !defined _COMPLEX_SYMMETRIC_
	     // complex-Hermitian case
	     GEMM("C","N",&n_size,&n_size,&m_size,
		  &alpha,
		  pvBLL,&m_size,
		  pvBLinvL,&m_size,
		  &beta,
		  pvBLinvD,&n_size,1,1);
#else
	     GEMM("T","N",&n_size,&n_size,&m_size,
		  &alpha,
		  pvBLL,&m_size,
		  pvBLinvL,&m_size,
		  &beta,
		  pvBLinvD,&n_size,1,1);
#endif
	  } // end if
	  
	  // successively downdate "n" by the size "n_size" of the diagonal block
	  sumn-=n_size;
	  
	  // symmetrize BLinv{k}.D
	  for (j=0; j<n_size; j++) {
	      // advance to the strict lower triangular part of column j
	      pvBLinvD+=j+1;

	      // pointer to BLinv{k}.D(j,j+1)
	      pv=pvBLinvD+n_size-1;
	      // symmetrize data from the strict lower/upper triangular part of row j
	      for (i=j+1; i<n_size; i++, pv+=n_size) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		  val=(*pv+*pvBLinvD)/2.0;
		  *pv=*pvBLinvD++=val;
#else
		  val.r=(pv->r+      pvBLinvD->r )/2.0;
		  val.i=(pv->i+CONJG(pvBLinvD->i))/2.0;
		  *pv=val;
		  pvBLinvD->r=val.r; pvBLinvD->i=CONJG(val.i);
		  pvBLinvD++;
#endif
	      } // end for i
	  } // end for j
#ifdef PRINT_INFO
	  printf("MYSYMSELBINV: %d-th inverse diagonal block computed\n", k+1);fflush(stdout);
	  pvBLinvD=BLinv.valD[k];
	  printf("        ");
	  for (j=0; j<n_size; j++)
	      printf("%8d", piBLinvJ[j]+1);
	  printf("\n");fflush(stdout);
	  for (i=0; i<n_size; i++) {
	      printf("%8d", piBLinvJ[i]+1);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvD[i+j*n_size]);
#else
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvD[i+j*n_size].r);
	      printf("\n");
	      for (j=0; j<n_size; j++)
		  printf("%8.1le", pvBLinvD[i+j*n_size].i);
#endif
	      printf("\n");
	      fflush(stdout);
	  }
#endif

          k--;
    } // end while k>=0

    // release auxiliary memory 
    free(gemm_buff);
    free(block);

    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // --------------------------   end main loop   ---------------------------
    //  now we have computed a selected inverse BLinv ~ BL^{-*} BiD BL^{-1}
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
#ifdef _PROFILING_
    time_selbinv=omp_get_wtime()-timeBegin;
#endif

	  




#ifdef _PROFILING_
    timeBegin = omp_get_wtime();
#endif
    iA->nc=iA->nr=n;
    iA->ncol  =(integer *) calloc((size_t)n,sizeof(integer));
    iA->rowind=(integer **)malloc((size_t)n*sizeof(integer *));
    iA->val   =(FLOAT **)  malloc((size_t)n*sizeof(FLOAT *));

    nz=0;
    for (k=0; k<nblocks; k++) {
        n_size=BLinv.nblockcol[k];
        m_size=BLinv.nblockrow[k];
	piBLinvJ=BLinv.colind[k];
	piBLinvI=BLinv.rowind[k];

	// allocate memory for lower triangular part of the permuted matrix
	for (r=0; r<n_size; r++) {
	    // column j
	    j=piBLinvJ[r];
	    q=perm[j];
	    s=m_size+n_size-r;
	    iA->ncol[q]=s;
	    iA->rowind[q]=(integer *)malloc((size_t)s*sizeof(integer));
	    iA->val[q]   =(FLOAT *)  malloc((size_t)s*sizeof(FLOAT));
	    // count number of nonzeros
	    nz+=s;
	} // end for j
	
	pvBLinvD=BLinv.valD[k];
	pvBLinvL=BLinv.valE[k];
	for (r=0; r<n_size; r++) {
	    // column j
	    j=piBLinvJ[r];
#ifdef PRINT_CHECK
	    if (j<0 || j>=n) {
	       printf("BLinv{%ld}.D, column index %ld out of range\n",k+1,j+1);
	       fflush(stdout);
	    }
#endif
	    q=perm[j];
	    piiAI=iA->rowind[q];
	    pviAv=iA->val[q];
	    // SL(q,q)
	    if (SL!=NULL)
	       deltaj=SL[q];
	    else
	       deltaj=1.0;
	    // column of BLinvD associated with j
	    pv=pvBLinvD+n_size*r;
	    // copy lower triangular part
	    l=0;
	    for (s=r; s<n_size; s++) {
	        i=piBLinvJ[s];
#ifdef PRINT_CHECK
		if (i<0||i>=n) {
		   printf("BLinv{%ld}.D, row index %ld out of range\n",k+1,j+1);
		   fflush(stdout);
		}
		if (perm[i]<0 || perm[i]>=n) {
		   printf("BLinv{%ld}.D, permuted row index %ld out of range\n",k+1,perm[i]+1);
		   fflush(stdout);
		}
		if (q<0 || q>=n) {
		   printf("BLinv{%ld}.D, column index %ld out of range\n",k+1,q+1);
		   fflush(stdout);
		}
#endif
		if (SL!=NULL)
		   deltai=SL[perm[i]];
		else
		   deltai=1.0;
	        piiAI[l]=perm[i];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		pviAv[l]=deltai*pv[s]*deltaj;
#else
		pviAv[l].r=deltai*pv[s].r*deltaj;
		pviAv[l].i=deltai*pv[s].i*deltaj;
#endif
		l++;
	    } // end for s
	    
	    // column of BLinvL associated with j
	    pv=pvBLinvL+m_size*r;
	    for (s=0; s<m_size; s++) {
	        i=piBLinvI[s];
#ifdef PRINT_CHECK
		if (i<0 || i>=n) {
		   printf("BLinv{%ld}.L, row index %ld out of range\n",k+1,i+1);
		   fflush(stdout);
		}
		if (perm[i]<0 || perm[i]>=n) {
		   printf("BLinv{%ld}.L, permuted row index %ld out of range\n",k+1,perm[i]+1);
		   fflush(stdout);
		}
		if (q<0 || q>=n) {
		   printf("BLinv{%ld}.L, column index %ld out of range\n",k+1,q+1);
		   fflush(stdout);
		}
#endif
		if (SL!=NULL)
		   deltai=SL[perm[i]];
		else
		   deltai=1.0;
	        piiAI[l]=perm[i];
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		pviAv[l]=deltai*pv[s]*deltaj;
#else
		pviAv[l].r=deltai*pv[s].r*deltaj;
		pviAv[l].i=deltai*pv[s].i*deltaj;
#endif
		l++;
	    } // end for s
	} // end for r

	// block column k is not longer required
	free(pvBLinvD);
	if (pvBLinvL!=NULL) free(pvBLinvL);
    } // end for k
    iA->nnz=nz;

    // BLinv is not longer required
    free(BLinv.valD);
    free(BLinv.valE);
#ifdef PRINT_INFO
    printf("MYSYMSELBINV: sparse list nz=%ld prepared\n",iA->nnz);fflush(stdout);
#endif

    

    // count number of nz for missing triangular part
    ncol=(integer *)calloc(n,sizeof(integer));
    for (j=0; j<n; j++) {
        // number of nonzeros in column j
        k=iA->ncol[j];
	// nonzero index list in column j
        piiAI=iA->rowind[j];
	for (l=0; l<k; l++) {
	    // row index i of iA(i,j)
            i=piiAI[l];
	    // additional entry in column i
	    if (i!=j)
	       ncol[i]++;
	} // end for l
    } // end for k

    // allocate additional memory
    for (j=0; j<n; j++) {
        // number of nonzeros in column j
        k=iA->ncol[j];
	// nonzero index list in column j
        piiAI=iA->rowind[j];
	// associated nonzeros
        pviAv=iA->val[j];
	// additional number of nonzeros
	q=k+ncol[j];
	iA->nnz+=ncol[j];
        piiAI=iA->rowind[j]=(integer *)realloc(piiAI,(size_t)q*sizeof(integer));
        pviAv=iA->val[j]   =(FLOAT *)  realloc(pviAv,(size_t)q*sizeof(FLOAT));
    } // end for k
#ifdef PRINT_INFO
    printf("MYSYMSELBINV: nz=%ld for output matrix defined\n",iA->nnz);fflush(stdout);
#endif
    // reset ncol
    for (j=0; j<n; j++)
        ncol[j]=0;
    // insert additional entries
    for (j=0; j<n; j++) {
        // number of nonzeros in column j
        k=iA->ncol[j];
	// nonzero index list in column j
        piiAI=iA->rowind[j];
	// associated nonzeros
        pviAv=iA->val[j];
	for (l=0; l<k; l++) {
	    // row index i of iA(i,j)
            i=piiAI[l];
	    // additional entry in column i
	    if (i!=j) {
	       q=iA->ncol[i]+ncol[i];
	       iA->rowind[i][q]=j;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       iA->val[i][q]=pviAv[l];
#else
	       iA->val[i][q].r=pviAv[l].r;
	       iA->val[i][q].i=CONJG(pviAv[l].i);
#endif
	       ncol[i]++;
	    } // end if
	} // end for l
    } // end for j
    // update number of nonzeros
    for (j=0; j<n; j++) {
        // number of nonzeros in column j
        k=iA->ncol[j];
	// additional number of nonzeros
	q=k+ncol[j];
	// new number of nonzeros
        iA->ncol[j]=q;
    } // end for k
    // sort columns
    for (j=0; j<n; j++) {
        // number of nonzeros in column j
        q=iA->ncol[j];
	// nonzero index list in column j
        piiAI=iA->rowind[j];
	// associated nonzeros
        pviAv=iA->val[j];
	QSORT(pviAv,piiAI,ncol,&q);
    } // end for j
#ifdef PRINT_INFO
    printf("MYSYMSELBINV: output matrix reordered\n");fflush(stdout);
#endif

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
    iA->isreal=1;
#else
    iA->isreal=0;
#endif
    
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_
    iA->issingle=1;
#else
    iA->issingle=0;
#endif
    
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
    iA->issymmetric=iA->ishermitian=0;
#else
#ifdef _COMPLEX_SYMMETRIC_
    iA->issymmetric=0; iA->ishermitian=0;
#else
    iA->issymmetric=0; iA->ishermitian=0;
#endif
#endif
    
    iA->isskew=0;
    
    free(ncol);
#ifdef _PROFILING_
    time_conversion=omp_get_wtime()-timeBegin;
#endif
       
#ifdef PRINT_INFO
    printf("MYSYMSELBINV: memory released\n");fflush(stdout);
#endif
    
#ifdef _PROFILING_
    time_symselbinv=omp_get_wtime()-time_symselbinv;
    printf("profiling summary\n");
    printf("1) Block selected inversion                 %12.4le\n",time_selbinv);
    printf("2) conversion block->scalar                 %12.4le\n",time_conversion);
    printf("3) remaining parts                          %12.4le\n",MAX(0.0,time_symselbinv-time_selbinv-time_conversion));
    printf("Total MYSYMSELBINV time %12.4le\n\n",time_symselbinv);

    fflush(stdout);
#endif
    
    return;
} // end symselbinv

