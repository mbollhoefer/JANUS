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

#define IABS(A)         (((A)>=(0))?(A):(-(A)))
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#define  SIZE_D    4

void PCG(integer *n, FLOAT  *rhs, FLOAT  *sol, integer *ipar, REALS *fpar,
	 FLOAT  *w)
{

  /* $Id: pcg.c 6094 2020-03-28 16:17:50Z bolle $
-----------------------------------------------------------------------
     This is a implementation of the PRECONDITIONED Conjugate Gradient
     (CG) method for solving linear system.

-----------------------------------------------------------------------

     n     size of the system

     rhs   right hand side

     sol   on entry, a given approximate solution
           on exit, the new approximate solution

     ipar  integer parameters
   
     fpar  real parameters
 
     w     work space

-----------------------------------------------------------------------

     ipar( 1):  status flag. depending on its values the algorithm is quit 
                for different requests
                0 initialization/return
                1 matrix-vector multiplication y=Ax
                2 transposed matrix-vector multiplication y=A^Tx
                3 left preconditioning
                4 left preconditioning transposed
                5 right preconditioning
                6 right preconditioning transposed
               10 test external stopping routine
               -1 too many iteration steps
               -2 not enough work space
               -3 algorithm breaks down
                  for Arnoldi-type methods refer to ipar(12) for a
                  more detailed description of the break-down

     ipar( 2):  0 no preconditioning used
                1 ONLY left preconditioning
                2 ONLY right preconditioning
                3 left AND right preconditioning

     ipar( 3):  stopping criteria
                  1  use residual norm as computed by cg
                  2  use rel. tol. * ||b|| + abs. tol.
                999  external stopping criterium, 

     ipar( 4):  size of working array w

     ipar( 5):  restart length (only used for GMRES,FOM,...

     ipar( 6):  maximum number of iteration steps

     ipar( 7):  counter for the number of steps (matrix-vector multiplications)

     ipar( 8):  address of the source x when applying one of the operations
                y in {Ax, A^Tx, Ml^{-1}x, Mr^{-1}x}

     ipar( 9):  associated target address for y

     ipar(10):  status flag. depending on its values the algorithm is 
                reentered at different labels
                0 initialization
                1 return after INITIAL matrix-vector multiplication y=Ax
                2 return after preconditioning
                3 return after matrix-vector multiplication y=Ax

     ipar(11):  flag for successful external stopping criterion
                 0 failure
                 1 success

     ipar(12):  for Arnoldi-type methods: break down, reason
		 -1  due to zero input vector
		 -2  since input vector contains abnormal numbers
		 -3  since input vector is a linear combination of others
		 -4  since triangular system in GMRES/FOM/etc. has null rank
                 
     ipar(13):
     ipar(14):
     ipar(15):
     ipar(16):


     fpar( 1):  relative error tolerance

     fpar( 2):  absolute error tolerance
     fpar( 3):  initial bound of the residual
     fpar( 4):  target bound for the residual
     fpar( 5):  current error bound of the residual
     fpar( 6):  current error bound of the delta x  (or the residual)
     fpar( 7):  convergence rate, here: norm of the preconditioned residual squared
     fpar( 8):
     fpar( 9):
     fpar(10):
     fpar(11): flop counter
     fpar(12):
     fpar(13):
     fpar(14):
     fpar(15):
     fpar(16):


     fpar(7) is used here internally to store <r, r>.
     w(:,1) -- residual vector
     w(:,2) -- P, the conjugate direction
     w(:,3) -- A P, matrix multiply the conjugate direction
     w(:,4) -- temporary storage for results of preconditioning
     w(:,5) -- change in the solution (sol) is stored here until
               termination of this solver

-----------------------------------------------------------------------
     Authors:

	Matthias Bollhoefer, TU Braunschweig

     Date:

       code modified from SPARSKIT-cg of Yousef Saad, adapted for PCG 
       and the complex case
       August 23, 2003. ILUPACK V1.0

     Notice:

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

     Availability:

	This file is located at

	http://ilupack.tu-bs.de/
-----------------------------------------------------------------------
  */

  // local variables

  static integer i, j, nupos;
  static REALS alpha, rho, nrmx, nrmx0, nrmd, nrm0, sumnu, nu[SIZE_D]; 
  static FLOAT falpha;
  static integer lp,rp;

  // check the status of the call

  if (ipar[0]<=0) ipar[9]=0;

  switch (ipar[9]) {
  case 1:
    goto label10;
    break;
  case 2: 
    goto label30;
    break;
  case 3: 
    goto label60;
    break;
  case 4: 
    goto label70;
    break;
  } // end switch

  // initialization
  i=5**n;
  j=1;
  BISINIT(ipar,fpar,&i,&j,&lp,&rp,w);
  if (ipar[0]<0) return;
  // position inside array nu
  nupos=1;
  nrmx0=-1.0;

  // ------------------------------------------------------------------
  // request for matrix vector multiplication y=A*x for initialization

  // set status flag for a request for performing y=A*x
  ipar[0]=1;
  // where does the source x start: w(:,2)
  ipar[7]=*n+1;
  // where do we ask for the return value y: w(:,3)
  ipar[8]=ipar[7]+*n;
  // return flag, indicate that we want to re-enter at label 10
  ipar[9]=1;
  // copy the initial solution to the buffer for the source x
  // w(:,2)=sol
  i=1;
  COPY(n,sol,&i,&w[0+*n],&i);
  return;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // we have performed a matrix-vector multiplication
  // the result y=Ax is returned in w(:,3)
  // now compute the initial residual r=b-Ax
     
  // increment step counter
  label10:
    ipar[6]=ipar[6]+1;
    ipar[12]=1;
    // preparations if the relative error in the energy norm is required
    if (ipar[2]==3) {
       i=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
       // scalar product (b,x)
       falpha=DISTDOT( n,rhs,&i,sol,&i);
       // real part
       alpha=2*falpha;
#else
       // scalar product (b,x)
#ifdef _USE_MKL_
       DISTDOT(&falpha,n,rhs,&i,sol,&i);
#else
       falpha=DISTDOT( n,rhs,&i,sol,&i);
#endif       
       // real part
       alpha=2*falpha.r;
#endif
       // 2* real(b,x_0) - ||x_0||_A^2
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
       nrm0=alpha-DISTDOT( n,sol,&i,&w[0+2**n],&i);
#else
#ifdef _USE_MKL_
       DISTDOT(&falpha,n,sol,&i,&w[0+2**n],&i);
       nrm0=alpha-falpha.r;
#else
       falpha=DISTDOT( n,sol,&i,&w[0+2**n],&i);
       nrm0=alpha-falpha.r;
#endif 
#endif 
       for (i=0; i<SIZE_D; i++)
	   nu[i]=0.0;
       // flop counter
       fpar[10]=fpar[10]+4**n;
    } // end if
    // build residual r=b-Ax, store it in w(:,1)
    for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
        w[i]  =rhs[i]  -w[i+2**n];
#else
        w[i].r=rhs[i].r-w[i+2**n].r;
        w[i].i=rhs[i].i-w[i+2**n].i;
#endif
    } // end for i
    // flop counter
    fpar[10]=fpar[10]+*n;
    // norm of the initial residual ||r||
    i=1;
    fpar[2]=NRM(n,w,&i);
    // flop counter
    fpar[10]=fpar[10]+2**n;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // current norm, here: initial residual
  fpar[4]=fpar[2];
  // which stopping criterion did we select
  // case 1: reltol * ||b|| + abstol
  if (IABS(ipar[2])==2) {
     // error bound fpar(1)=relative tolerance
     //             fpar(2)=absolute tolerance
     // model: reltol*||b||+abstol
     i=1;
     fpar[3]=fpar[0]*sqrt(DISTDOT2(n,rhs,&i))+fpar[1];
     // increase flop counter
     fpar[10]=fpar[10]+2**n;
  } // end if
  //case 2: reltol * ||r|| + abstol
  else if (ipar[2]!=999) {
     // error bound fpar(1)=relative tolerance
     //             fpar(2)=absolute tolerance
     // model: reltol*||r||+abstol, where r is the residual
     fpar[3]=fpar[0]*fpar[2]+fpar[1];
     //fpar(4) not used if ipar(3)=3
  } // end if-else if
  // ------------------------------------------------------------------



  // ------------------------------------------------------------------
  // ---------------------     START MAIN LOOP     --------------------
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // if preconditioning z=M^{-1}r is desired 

  label20:
    if (lp || rp) {
       // set status flag for a request for preconditioning
       if (lp) 
	  ipar[0]=3;
       else
	  ipar[0]=5;
       // where does the source r start: w(:,1)
       ipar[7]=1;
       // target z for preconditioning z=M{-1}r: w(:,4)
       ipar[8]=3**n+1;
       // return flag, indicate to return to label 30 when re-entering 
       // the subroutine
       ipar[9]=2;
       return;
    } // end if
    else {
       // use w(:,4)=w(:,1)
       i=1;
       COPY(n,w,&i,&w[0+3**n],&i);
    } // end if-else
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // compute scalar product rho=r^Hz of between the residual r and the 
  // preconditioned residual z. rho must be real and positive
  label30:
    i=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    rho=   DISTDOT( n,w,&i,&w[0+3**n],&i);
#else
#ifdef _USE_MKL_
    DISTDOT(&falpha,n,w,&i,&w[0+3**n],&i);
    rho=falpha.r;
#else
    falpha=DISTDOT( n,w,&i,&w[0+3**n],&i);
    rho=falpha.r;
#endif 
#endif 
    
    // increment flop counter
    fpar[10]=fpar[10]+2**n;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // init/update conjugate search direction
  if (ipar[6]==1) {
     // initial conjugate direction stored in w(:,2)
     // p=z
     i=1;
     COPY(n,&w[0+3**n],&i,&w[0+*n],&i);
  } // end if
  else {
    // parameter for update of the search direction
    // r_new^Hz_new / r_old^Hz_old
    alpha=rho/fpar[6];
    // p = z + alpha p
    for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
        w[i+*n]  =w[i+3**n]  +alpha*w[i+*n];
#else
        w[i+*n].r=w[i+3**n].r+alpha*w[i+*n].r;
        w[i+*n].i=w[i+3**n].i+alpha*w[i+*n].i;
#endif
    } // end for i
    // increment flop counter
    fpar[10]=fpar[10]+2**n;
  } // end if-else
  // use fpar(7) to store the actual rho
  fpar[6]=rho;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // before iteration can continue, we need to compute A * p, where
  // p=w(:,2) is the conjugate search direction
  // set status flag as request for matrix-vector multiplication w=Ap
  label40:
    ipar[0]=1;
    // source address: w(:,2), conjugate search direction p
    ipar[7]=*n+1;
    // address of the target w: w(:,3)
    ipar[8]=2**n+1;
    // return flag, indicate to return to label 60 when re-entering the
    // subroutine
    ipar[9]=3;
    return;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // continueing with the iterations. 

  // advance counter for the matrix-vector multiplications
  label60:
    ipar[6]=ipar[6]+1;
    // compute alpha=p^HAp
    i=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha=DISTDOT(  n,&w[0+*n],&i,&w[0+2**n],&i);
#else
#ifdef _USE_MKL_
    DISTDOT(&falpha,n,&w[0+*n],&i,&w[0+2**n],&i);
    alpha=falpha.r;
#else
    falpha=DISTDOT( n,&w[0+*n],&i,&w[0+2**n],&i);
    alpha=falpha.r;
#endif 
#endif 
    // increment flop counter
    fpar[10]=fpar[10]+2**n;
    // alpha=0 or +/-oo ?
    if (BRKDN(&alpha,ipar)) goto label900;
    // parameter to update approx. solution and residual
    // dx = dx + (r^Hz/p^HAp)  p, dx: sol+w(:,5) approx. sol.
    // r  = r  - (r^Hz/p^HAp) Ap,  r: w(:,1)     residual
    //                             p: w(:,2)     conjg. search direction
    //                            Ap: w(:,3)     A * conjg. search direction
    // dx = dx + alpha  p
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    falpha  =rho/alpha;
#else
    falpha.r=rho/alpha;
    falpha.i=0.0;
#endif
    if (ipar[2]==3) {
       // additional contribution to the energy norm estimate
       nu[nupos-1]=rho*rho/alpha;
       // compute nu_j,d
       sumnu=0.0;
       for (i=0; i<SIZE_D; i++)
	   sumnu=sumnu+nu[i];
       fpar[4]=sqrt(fabs(sumnu));
       // update nu_0,j + 2* real(b,x)-||x_0||_A^2
       nrm0=nrm0+nu[nupos-1];
       nupos=nupos+1;
       if (nupos>SIZE_D) nupos=1;
    } // end if
    i=1;
    AXPY(n,&falpha,&w[0+*n],&i,&w[0+4**n],&i);
    // r = r   - alpha Ap
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    falpha  =-falpha;
#else
    falpha.r=-falpha.r;
    falpha.i=-falpha.i;
#endif
    AXPY(n,&falpha,&w[0+2**n],&i,w,&i);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha=-falpha;
#else
    alpha=-falpha.r;
#endif
    // increment flop counter
    fpar[10]=fpar[10]+4**n;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // are we ready to terminate ?

  // external stopping criterion
  if (ipar[2]==999) {
     // set status flag for a request for external stopping criterion
     ipar[0]=10;
     // first address w(:,5) (approx. solution dx)
     ipar[7]=4**n+1;
     // second address w(:,4), some target
     ipar[8]=3**n+1;
     // return flag, indicate to return to label 70 when re-entering the
     // subroutine
     ipar[9]=4;
     return;
  } // end if
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // check stopping criteria
  label70:
    if (ipar[2]==999) {
       if (ipar[10]==1) goto label900;
       // w(:,1) residual
       // w(:,2) conjg. search direction
       // alpha*w(:,2) ~ delta x change in the solution
    } // end if
    else if (ipar[2]==3 && ipar[6]>SIZE_D) {
       if ((ipar[6]>SIZE_D && sumnu<fpar[0]*fpar[0]*(nrm0+sumnu))
	   || ipar[6]>=ipar[5] || alpha==0.0) {
	  fpar[3]=fpar[0]*sqrt(nrm0+sumnu);
	  // at first glance it looks as if we have succeeded
	  // but to be save we check that the update has a small
	  // relative norm compared with x
	  // ||p|| <= reltol ||x||/2
	  // for ||x|| we use ||x0||+||dx|| as estimate
	  i=1;
	  nrmd=NRM(n,&w[0+*n],  &i);
	  nrmx=NRM(n,&w[0+4**n],&i);
	  if (nrmx0<0.0) nrmx0=NRM(n,sol,&i);
	  // increase flop counter
	  ipar[10]=ipar[10]+4**n;
	  if ((ipar[6]>=ipar[5]) || nrmd<=0.5*fpar[0]*(nrmx+nrmx0)) goto label900;
       } // end if
    } // end else if
    else if (ipar[2]!=3 || (ipar[2]==3 && ipar[6]<=SIZE_D)) {
       i=1;
       if (STOPBIS(n,ipar,&i,fpar,w,&w[0+*n],&alpha)) goto label900;
    } // end if-else if-else if
  // ------------------------------------------------------------------


  // continue the iterations
  goto label20;
  // ------------------------------------------------------------------
  // ----------------------     END MAIN LOOP     ---------------------
  // ------------------------------------------------------------------



  // ------------------------------------------------------------------
  // clean up 

  // sol = sol + w(:,5)
  label900:
    TIDYCG(n,ipar,fpar,sol,&w[0+4**n]);
    return;
} // end-of-pcg
