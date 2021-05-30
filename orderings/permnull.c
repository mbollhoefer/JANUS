#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <janus.h>

#include <ilupackmacros.h>
#define MAX_LINE        255
#define STDERR          stderr
#define STDOUT          stdout
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))





/*
 scaling and permutation driver, here dummy driver that does at most scaling

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


 In particular here no permutation is computed. nB=n
 */
integer PERMNULL(CSRMAT A, FLOAT *prowscale, FLOAT *pcolscale, 
	     integer *p, integer *invq, integer *nB,  ILUPACKPARAM *param)

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


   interface driver written by Matthias Bollhoefer, July/November 2003/Feb 2014
 */

   integer i,j,k,scale,nrm,ierr=0;


   // start by rescaling the system, param defines the norms
#include "scaleprefix.c"

   // build both permutations
   for (i=0; i<A.nc; i++) {
       p[i]=invq[i]=i+1;
   } // end for

   return (ierr);
} /* end permnull */
