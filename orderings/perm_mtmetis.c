#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include <mtmetis.h>
#include <janus.h>

#include <ilupackmacros.h>

#define MAX_LINE        255
#define STDERR          stderr
#define STDOUT          stdout
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))

// #define PRINT_INFO



/*
 scaling and permutation driver, here for multilevel (edge) nested dissection
 with minimum degree from METIS (Karypis & Kumar)

 Given an n x n matrix A in compressed sparse row format this routine computes
 row and column scalings as well as row and column permutations such that

 A -> Dr * A * Dc,
 
 where Dr and Dc are diagonal matrices stored in the vectors proscal and 
 pcoscal. The scaling is explicitly applied to A.

 The permutation p,invq refer to permutation matrices P and Q^{-1} such  that

  P^T*(Dr * A * Dc)*Q  

 is hopefully better suited for the ILU.

 The routine returns a leading block size nB<=n
                        /B F\ 
  P^T*(Dr * A * Dc)*Q = |   | 
                        \E C/
 In this case the permutation recommends that the ILU is only applied to the
 leading block B of size nB x nB


 In particular MMD uses a symmetric reordering that reorders the system
 such that nodes from the elimination graph with lowest degree precede
 those with higher degree which may be advantageous for direct solvers. nB=n
 */
integer PERMMTMETIS(CSRMAT A, FLOAT *prowscale, FLOAT *pcolscale,
		    integer *p, integer *invq, integer *nB, ILUPACKPARAM *param)

{
/*
    A          n x n matrix in compressed sparse row format
               A is altered by the routine since the scalings are directly
               applied to A. Hopefully the scaling are only powers of 2,
               which is the case for all scaling routines from ILUPACK.

    proscal,   vectors of length n that store the row and column scalings Dr
    pcoscal    and Dc
               A -> Dr*A*Dc

    p,invq     permutation vectors p and q^{-1} of length n which refer to
               row / column permutation matrices P and Q^{-1}.
               P^T*(Dr * A * Dc)*Q hopefully better suited for ILU

    nB         leading blocksize nB
                                     /B F\ 
               P^T*(Dr * A * Dc)*Q = |   |,  B matrix of size nB x nB
                                     \E C/
               nB is the recommended block size for the application of an ILU

    param      ILUPACK parameters
               param->ipar[7] gives information about the requested scaling

               param->ipar[7] &  512       indicates initial preprocessing,
                                           initial permutation routine in
                                           the main ILU driver called AMGFACTOR

               param->ipar[7] & 1024       indicates regular reordering,
                                           the second permutation routine in
                                           the main ILU driver called AMGFACTOR

               param->ipar[7] &  512+1024  indicates final pivoting,
                                           the third permutation routine in
                                           the main ILU driver called AMGFACTOR


               depending on the circumstances (initial preprocessing, regular
               reordering, final pivoting)

               param->ipar[7] & ( 1+  2+  4)   initial preprocessing

               param->ipar[7] & ( 8+ 16+ 32)   regular reordering

               param->ipar[7] & (64+128+256)   final pivoting

               give information on whether row scaling (lowest bit 1,8,64), 
               column scaling (medium bit 2,16,128) should be applied.
               The highest bit (4,32,256) defines the order in which the 
               scalings should be perfomed. If the highest bit is cleared, then
               we start with row scaling.

    
               param->ipar[8] defines the norm that should be used

               param->ipar[8] & (   1+   2+   4+   8+   16)  initial preprocessing

               param->ipar[8] & (  32+  64+ 128+ 256+  512)  regular reordering

               param->ipar[8] & (1024+2048+4096+8192+16384)  final pivoting

               The five bits (values [0,...,31] up to shifts) are used for
               0,1,2          infinity norm, 1-norm, 2-norm
               3              spd scaling using the square root of the diagonal
                              entries

               The scaling routines that are defined with ILUPACK only use
               the nearest powers of 2 for scaling


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        July/November 2003, 2005. ILUPACK V2.1

    Notice:

	Copyright (c) 2005 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
 */

   integer i,j,k,l,scale,nrm,ierr=0, numflag, options[10], iflag, *vwgt;
   double *mtmetis_options=NULL;
   size_t mem;
   CSRMAT B;
   

   // start by rescaling the system, param defines the norms
#include "scaleprefix.c"


   // allocate memory with respect to the specific permutation routine
   param->nibuff=MAX(param->nibuff,(size_t)A.nc+1);
   param->ibuff=(integer *)   REALLOC(param->ibuff,param->nibuff*sizeof(integer),
				      "perm_mtmetis:ibuff");


   // build B=|A|+|A^T|
   // param->ibuff will also be used for B.ia!
   // remember that the leading 2n are already in use for symrcm
   // param->iaux will be used for B.ja if possible
   SETUPGRAPH(A,&B, param->ibuff, param->iaux,param->niaux);

   // check if MALLOC was needed to set up B.ja
   if (B.ia[B.nc]<0) {
      iflag=0;
      B.ia[B.nc]=-(B.ia[B.nc]);
   }
   else  {// nz+B.ia[B.nc]-1 entries from iaux were already in use 
      iflag=B.ia[B.nc]-1; 
   }

   // reorder the system using multilevel nested dissection & minimum degree
   // get rid of the diagonal entries
   k=0;
   for (i=0; i<A.nc; i++) {
       l=k;
       for (j=B.ia[i]-1; j<B.ia[i+1]-1;j++) {
	   // diagonal entry
	   if (B.ja[j]==i+1) {
	      // skip diagonal entry
	      k++;
	   }
	   else {
	      // shift diagonal entry
	      B.ja[j-k]=B.ja[j];
	   } // end if
       } // end for j
       B.ia[i]-=l;
   } // end for i
   B.ia[i]-=k;


   // reorder the system using multilevel nested dissection & minimum degree
   // switch from FORTRAN style to C-style
   for (i=0; i<=A.nc; i++)
       B.ia[i]--;
   for (i=0; i<B.ia[A.nc]; i++)
       B.ja[i]--;
   // default options
#ifndef _FULL_MTMETIS_
   int mt=omp_get_max_threads();
   omp_set_num_threads(MY_MTMETIS_THREADS);
#endif
   mtmetis_options=mtmetis_init_options();

   // matrix is non-diagonal
   if (B.ia[A.nc]>1 && A.nc>16) {
      vwgt=(integer *)malloc((size_t)A.nc*sizeof(integer));
      for (i=0; i<A.nc; i++)
	  vwgt[i]=1;
#ifdef PRINT_INFO
      printf("Call MT-Metis\n");fflush(stdout);
      for (i=0; i<A.nc; i++) {
	  printf("column %ld:\n",i);
	  for (j=B.ia[i]; j<B.ia[i+1]; j++)
	      printf("%4ld",B.ja[j]);
	  printf("\n");
	  fflush(stdout);
      } // end for i
      printf("\np:");
      for (i=0; i<A.nc; i++)
	  printf("%4ld",p[i]);
      printf("\ninvq:");
      for (i=0; i<A.nc; i++)
	  printf("%4ld",invq[i]);
      printf("\n");
      fflush(stdout);
#endif
      MTMETIS_NodeND(&A.nc,B.ia,B.ja, vwgt,mtmetis_options,invq,p);
#ifdef PRINT_INFO
      printf("MT-Metis completed\n");fflush(stdout);
#endif
      for (i=0; i<A.nc; i++) {
	  p[i]++;
	  invq[i]++;
      } // end for i
      free(vwgt);
   }
   else
      for (i=0; i<A.nc; i++) 
	  p[i]=invq[i]=i+1;

   // switch from C-style to FORTRAN style 
   for (i=0; i<B.ia[A.nc]; i++)
       B.ja[i]++;
   for (i=0; i<=A.nc; i++)
       B.ia[i]++;
   if (!iflag)
      B.ja=FREE(B.ja);
   free(mtmetis_options);
#ifndef _FULL_MTMETIS_
   omp_set_num_threads(mt);
#endif

   *nB=A.nc;


   return (ierr);
}
