#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <janus.h>

#include <ilupackmacros.h>

#define MAX_LINE        255
#define STDERR          stdout
#define STDOUT          stdout
#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))






/*
 scaling and permutation driver, here for indPP, a symmetric variant of
 indPQ (Yousef Saad)

 Given an n x n matrix A in compressed sparse row format this routine computes
 row and column scalings as well as row and column permutations such that

 A -> Dc * A * Dc,
 
 where Dc is a diagonal matrix stored in the vectors proscal and pcoscal. 
 The scaling is explicitly applied to A.

 The permutation p,invq refer to permutation matrices P and P^{-1} such  that

  P^T*(Dc * A * Dc)*P

 is hopefully better suited for the ILU.

 The routine returns a leading block size nB<=n
                        / B  F\ 
  P^T*(Dc * A * Dc)*P = |     | 
                        \F^T C/
 In this case the permutation recommends that the ILU is only applied to the
 leading block B of size nB x nB


 In particular indPP tries to construct a leading diagonal block B that
 follows a compromise between sparsity and diagonal dominance, typically nB<n
 */
integer SPDPERMPP(CSRMAT A, FLOAT *prowscale, FLOAT *pcolscale, 
	   integer *p, integer *invq, integer *nB, ILUPACKPARAM *param)
{
/*
    A          n x n matrix in compressed sparse row format. Only the upper
               triangular part (up to a symmetric permutation) is stored.
               A is altered by the routine since the scalings are directly
               applied to A. Hopefully the scaling are only powers of 2,
               which is the case for all scaling routines from ILUPACK.

    proscal,   vectors of length n that store the row and column scalings Dr
    pcoscal    and Dc. Here only pcoscal is referenced, since the matrix is
               symmetric. To have a consistent interface, proscal is formally
               needed
               A -> Dc*A*Dc

    p,invq     permutation vectors p and p^{-1} of length n which refer to
               row / column permutation matrices P and P^{-1}.
               Here p(invq)=id is required since the problem is symmetric.
               Both permutation vectors are needed (for technical reasons)
               P^T*(Dc * A * Dc)*P hopefully better suited for ILU

    nB         leading blocksize nB
                                     / B  F\ 
               P^T*(Dc * A * Dc)*P = |     |,  B matrix of size nB x nB
                                     \F^T C/
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


   interface driver written by Matthias Bollhoefer, July/November 2003
 */

   integer i,j,k,scale,nrm,ierr=0;


   // start by rescaling the system, param defines the norms
   // the ILUPACK parameters have to be set up such that a symmetric rescaling
   // is performed
#include "scaleprefix.c"



   // allocate memory with respect to the specific permutation routine
   param->nibuff=MAX(param->nibuff,2*(size_t)A.nc);
   param->ibuff=(integer *)   REALLOC(param->ibuff,param->nibuff*sizeof(integer),
				  "permPP:ibuff");
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
   param->ndbuff=MAX(param->ndbuff,3*(size_t)A.nc);
   param->dbuff=(FLOAT *) REALLOC(param->dbuff,param->ndbuff*sizeof(FLOAT),
				  "permPP:dbuff");
#else
   // 6*n REALS <= 3*n COMPLEX spaces required 
   param->ndbuff=MAX(param->ndbuff,MAX(3,(6*sizeof(REALS))/sizeof(FLOAT))*(size_t)A.nc);
   param->dbuff=(FLOAT *) REALLOC(param->dbuff,param->ndbuff*sizeof(FLOAT),
				   "permPP:dbuff");
#endif

   // indPP permutation routine
   PPPERM(A, 50, invq, nB, 0.5,param->dbuff,param->ibuff);
   for (i=0; i<A.nc; i++) {
       p[invq[i]]=i;
   }
   for (i=0; i<A.nc; i++) {
       invq[i]++;
       p[i]++;
   }

   return (ierr);
}
