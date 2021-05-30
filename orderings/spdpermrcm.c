#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <sparspak.h>
#include <janus.h>

#include <ilupackmacros.h>

#define MAX_LINE        255
#define STDERR          stderr
#define STDOUT          stdout
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))





/*
 scaling and permutation driver, here for reverse Cuthill-McKee from SPARSPAK 
 (Alan George & Joseph Liu)

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


 In particular reverse Cuthill-McKee uses a symmetric reordering that enforces
 a variable band structure, in the sense of its graph, levels of neighbouring
 nodes are determined. nB=n
 */
integer SPDPERMRCM(CSRMAT A, FLOAT *prowscale, FLOAT *pcolscale,
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


   interface driver written by Matthias Bollhoefer, July/November 2003, 2005
 */

  integer i,j,k,bnd,bndt,scale,nrm,ierr=0, iflag;
   REALS sm,sm2, sigma, sigmat, x;
   CSRMAT B;


   // start by rescaling the system, param defines the norms
   // the ILUPACK parameters have to be set up such that a symmetric rescaling
   // is performed
#include "scaleprefix.c"


   // allocate memory with respect to the specific permutation routine
   param->nibuff=MAX(param->nibuff,3*(size_t)A.nc+1);
   param->ibuff=(integer *)   REALLOC(param->ibuff,param->nibuff*sizeof(integer),
				  "permrcm:ibuff");

   // param->ibuff+2n will also be used for B.ia!
   // remember that the leading 2n are already in use for symrcm
   // param->iaux will be used for B.ja if possible
   SETUPGRAPH(A,&B, param->ibuff+2*A.nc, param->iaux,param->niaux);

   // check if MALLOC was needed to set up B.ja
   if (B.ia[B.nc]<0) {
      iflag=0;
      B.ia[B.nc]=-(B.ia[B.nc]);
   }
   else  {// nz+B.ia[B.nc]-1 entries from iaux were already in use 
      iflag=B.ia[B.nc]-1; 
   }

   // call reverse Cuthill-McKee
   i=0;
   genrcm(&A.nc,B.ia,B.ja, p,param->ibuff,param->ibuff+A.nc); 
   if (!iflag)
      B.ja=FREE(B.ja);


   /* determine rows and columns that have extremely many nonzeros */
   for (i=0; i<A.nc; i++) 
       param->ibuff[i]=0;
   for (i=0; i<A.nc; i++) {
       for (j=A.ia[i]; j<A.ia[i+1]; j++) {
           /* column index A(i,k) */
           k=A.ja[j-1]-1;
	   // increment number of nonzeros except for the diagonal part 
	   if (k!=i) param->ibuff[k]++;
       }
   }
   // remember that only the upper triangular part is stored
   for (i=0; i<A.nc; i++) 
       param->ibuff[i]+=A.ia[i+1]-A.ia[i];

   /* arithmetic mean: nnz/n */
   x=2*(A.ia[A.nc]-1)/((REALS)A.nc);
   /* compute standard deviation */
   sm=0; sm2=0;
   j=0;
   for (i=0; i<A.nc; i++) {
       k=param->ibuff[i];
       sm+=k;
       sm2+=k*k;
   }
   /* standard deviation by row/column */
   sigma=sqrt((sm2+A.nc*x*x-2.0*x*sm)/A.nc); 
   
   /* now push those rows/columns to the end that have too many nonzeros */
   j=0;
   k=0;
   bnd =MAX(2*x,x+2*(sigma +.5));
   for (i=0; i<A.nc; i++) {
     if (param->ibuff[p[i]-1]>bnd) {
        invq[k]=p[i];
        k++;
     }
     else {
       p[j]=p[i];
       j++;
     }
   }
   /* reinsert skipped entries at the end */
   i=0;
   for (; j<A.nc; j++,i++)
       p[j]=invq[i];
   
   /* compute inverse permutation */
   for (i=0; i<A.nc; i++)
       invq[p[i]-1]=i+1;


   *nB=A.nc-k; 


   return (ierr);
} /* end spdpermrcm */
 
