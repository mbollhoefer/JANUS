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

#define IABS(A)          (((A)>=(0))?(A):(-(A)))
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)

#ifdef _DOUBLE_REAL_
#define MYSQMR       DSYMqmr
#define MYDISTDOT    ddot
#else
#define MYSQMR       SSYMqmr
#define MYDISTDOT    sdot
#endif



#else



#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)

#ifdef _SINGLE_COMPLEX_
#define MYSQMR       CSYMqmr
#define MYDISTDOT    cdotu
#else
#define MYSQMR       ZSYMqmr
#define MYDISTDOT    zdotu
#endif

#else

#ifdef _SINGLE_COMPLEX_
#define CONJG(A)     (conjg(A))
#define MYSQMR       CHERqmr
#define MYDISTDOT    cdotc
#else
#define CONJG(A)     (dconjg(A))
#define MYSQMR       ZHERqmr
#define MYDISTDOT    zdotc
#endif

#endif



#endif

void MYSQMR(integer *n,FLOAT *rhs,FLOAT *sol,integer *ipar,REALS *fpar,FLOAT *w)
{

  /* $Id: sqmr.c 6094 2020-03-28 16:17:50Z bolle $
    -----------------------------------------------------------------------
     SQMR: symmetric QMR method. Programmed with reverse
     communication, see the header for detailed specifications
     of the protocol.

     ipar(1) control flag for external operation
              1  y=Ax        matrix-vector multiplication

              3  z=ML^{-1} r left preconditioning

              5  z=MR^{-1} r right preconditioning

             10  external stopping criterion
     ipar(3) stopping criterion
              1 ||b-Ax_k|| <= reltol * ||b-Ax_0|| +abstol 
              2 ||b-Ax_k|| <= reltol * ||b|| +abstol 
              3 ||b-Ax_k|| <= reltol * (||b||+||A|| ||x_k||)
     in this routine, before successful return, the fpar's are
     fpar(3) == initial residual norm
     fpar(4) == target residual norm
     fpar(5) == current residual norm
     fpar(7) == current |rho| (rhok = <r, M^{-1}r>)
     fpar(8) == previous |rho| (rhokm1)

     w(:,1) -- r, the residual
     w(:,2) -- p, the projection direction
     w(:,3) -- y, a scratch vector to store A*p, or A*q.
     w(:,4) -- z, a scratch vector to store intermediate results
     w(:,5) -- x, updated solution
     w(:,6) -- d, update direction for x
     
     written by Matthias Bollhoefer, December 2004
     -----------------------------------------------------------------------
  */

  // local variables
  static integer delay, delaycount;
  static REALS theta, theta_old, c, c2, tau, adjust;
  static integer i,j;
  static FLOAT alpha, rho, rho_old, beta;
  static REALS aalpha, nrmx, nrmb, nrmd;
  static integer rp, lp;

  // status of the program

  if (ipar[0]<=0) ipar[9]=0;
  
  switch (ipar[9]) {
  case 1:
    goto label10;
    break;
  case 2: 
    goto label20;
    break;
  case 3: 
    goto label30;
    break;
  case 4: 
    goto label40;
    break;
  case 5: 
    goto label50;
    break;
  case 6: 
    goto label60;
    break;
  case 7: 
    goto label70;
    break;
  case 8: 
    goto label80;
    break;
  } // end switch

  // initialization, initial residual

  i=6**n;
  j=1;
  BISINIT(ipar,fpar,&i,&j,&lp,&rp,w);

  // adjust residual and preconditioned residual
  adjust=1.0;
  if (ipar[0]<0) return;
  
  // ------------------------------------------------------------------
  // request for matrix vector multiplication y = A x for initialization

  // set status flag for a request for performing y = A x
  ipar[0]=1;
  // where does the source x start: w(:,5)
  ipar[7]=4**n+1;
  // where do we ask for the return value y: w(:,3)
  ipar[8]=2**n+1;
  // return flag, indicate that we want to re-enter at first label 10
  ipar[9]=1;
  // copy the initial solution to the buffer for the source x=w(:,5)
  i=1;
  j=1;
  COPY(n,sol,&i,&w[0+4**n],&j);
  return;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // we have performed a matrix-vector multiplication
  // the result y=Ax is returned in w(:,3)
  // now compute the initial residual r=b-Ax
 
  // increment step counter
  label10:
    ipar[6]=ipar[6]+1;
    // write (6, '(A)') 'label 10'
    ipar[12]=ipar[12]+1;
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

    // left preconditioning
    if (lp) {
       // set status flag for a request for y = ML^{-1} r
       ipar[0]=3;
       // where does the source r start: w(:,1)
       ipar[7]=1;
       // where do we ask for the return value y: w(:,3)
       ipar[8]=2**n+1;
       //return flag, indicate that we want to re-enter at the second label 20
       ipar[9]=2;
       return;
    } // end if lp
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // we have computed the preconditioned residual y = ML^{-1} r
  // r is stored in w(:,1)
  // y = ML^{-1} r is stored in w(:,3)


  // set up parameters for QMR search directions
  label20:
    theta=0;
    // write (6, '(A)') 'label 20'

    // tau=||y||=||ML^{-1}r||
    i=1;
    if (lp) 
       tau=NRM(n,&w[0+2**n],&i);
    else
       tau=NRM(n,w,&i);


    // Now we set upt the search directions
    // p = MR^{-1} ML^{-1} r
    // request for right preconditioning p = MR^{-1} y
    if (rp) {
       // set status flag for a request for p = MR^{-1} y
       ipar[0]=5;
       // depending on whether we have left preconditioning in advance
       // the source vector is located at different places
       if (lp)
	  // where does the source y=ML^{-1}r start: w(:,3)
	  ipar[7]=2**n+1;
       else
	  // no left preconditioning, y = r = w(:,1)
	  ipar[7]=1;

       // where do we ask for the return value p: w(:,2)
       ipar[8]=*n+1;
       // return flag, indicate that we want to re-enter at the third label 30
       ipar[9]=3;
       return;
    } // end if 
    else { // !rp
       // p = y = ML^{-1} r
       i=1;
       j=1;
       // depending on whether we have left preconditioning in advance
       // the source vector is located at different places
       if (lp)
	  // p = y = ML^{-1} r
	  COPY(n,&w[0+2**n],&i,&w[0+*n],&j);
       else
	  // p = r
	  COPY(n,w,         &i,&w[0+*n],&j);
    } // end if-else rp
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // now we have computed the search direction p = MR^{-1} ML^{-1} r

  // r=w(:,1)
  // p=w(:,2)

  // rho = p^* r
  label30:
    i=1;
    j=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    rho=MYDISTDOT( n,&w[0+*n],&i,w,&j);
#else
#ifdef _USE_MKL_
    MYDISTDOT(&rho,n,&w[0+*n],&i,w,&j);
#else
    rho=MYDISTDOT( n,&w[0+*n],&i,w,&j);
#endif 
#endif 
    fpar[6]=FABS(rho);
    // flop counter
    fpar[10]=fpar[10]+2**n;
    // initial residual norm
    // write (6, '(A)') 'label 31'
    fpar[2]=NRM(n,w,&j);

    // current precnd. quasi residual norm
    fpar[4]=tau;
    // adjust different norms, for safety use additional factor 0.5
    if (fpar[2]!=0.0) adjust=0.5*(fpar[4]/fpar[2]);
    fpar[7]=1.0;
    j=1;
    // which stopping criterion did we select
    // case 2: reltol * ||b|| + abstol
    if (IABS(ipar[2])==2) {
       // error bound fpar(1)=relative tolerance
       //             fpar(2)=absolute tolerance
       // model: reltol*||b||+abstol
       // write (6, '(A)') 'label 32'
       nrmb=NRM(n,rhs,&j);
       fpar[3]=fpar[0]*nrmb+fpar[1];
       // increase flop counter
       fpar[10]=fpar[10]+2**n;
    } // end if
    // case 1: reltol * ||r|| + abstol
    else if (ipar[2]!=999 && ipar[2]!=3) {
       // error bound fpar(1)=relative tolerance
       //             fpar(2)=absolute tolerance
       // model: reltol*||r||+abstol, where r is the residual
       fpar[3]=fpar[0]*fpar[2]+fpar[1];
    } // end if-else if
    else if (ipar[2]==3) {
       // backward error
       // fpar(12): ||A||
       i=1;
       // write (6, '(A)') 'label 33'
       nrmb=NRM(n,rhs,&i);
       // write (6, '(A)') 'label 34'
       nrmx=NRM(n,sol,&i);
       delaycount=1;
       fpar[3]=adjust*fpar[0]*(nrmb+fpar[11]*nrmx);
    } // end if-else if-else if

    // did the method already converge
    if (ipar[2]>=0 && fpar[4]<=fpar[3]) {
       fpar[5]=fpar[2];
       goto label900;
    } // end if
    // write (6, '(A)') 'label 35'
  // ------------------------------------------------------------------



  // ------------------------------------------------------------------
  // ---------------------     START MAIN LOOP     --------------------
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // request for matrix vector multiplication y = A p

  // set status flag for a request for performing y = A p
  label300:
    ipar[0]=1;
    // where does the source p start: w(:,2)
    ipar[7]=*n+1;
    // where do we ask for the return value y: w(:,3)
    ipar[8]=2**n+1;
    // return flag, indicate that we want to re-enter at the fourth label 40
    ipar[9]=4;
    return;
  // ------------------------------------------------------------------


  // increment step counter
  label40:
    ipar[6]=ipar[6]+1;
    // write (6, '(A)') 'label 40'

    // alpha = p^* Ap
    i=1;
    j=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha=MYDISTDOT( n,&w[0+*n],&i,&w[0+2**n],&j);
#else
#ifdef _USE_MKL_
    MYDISTDOT(&alpha,n,&w[0+*n],&i,&w[0+2**n],&j);
#else
    alpha=MYDISTDOT( n,&w[0+*n],&i,&w[0+2**n],&j);
#endif
#endif
    // increase flop counter
    fpar[10]=fpar[10]+2**n;
    // check for break down
    aalpha=FABS(alpha);
    // |alpha|=0 or +/-oo ?
    if (BRKDN(&aalpha,ipar)) goto label900;

    // alpha = rho / p^* (A p)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha=rho/alpha;
#else
    aalpha=1.0/FABS2(alpha);
    beta.r= alpha.r*aalpha;
    beta.i=-alpha.i*aalpha;
    alpha.r=rho.r*beta.r-rho.i*beta.i;
    alpha.i=rho.r*beta.i+rho.i*beta.r;
#endif
    i=1;
    j=1;
    // r = r - alpha  Ap
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha  =-alpha;
#else
    alpha.r=-alpha.r;
    alpha.i=-alpha.i;
#endif
    AXPY(n,&alpha,&w[0+2**n],&i,w,&j);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha  =-alpha;
#else
    alpha.r=-alpha.r;
    alpha.i=-alpha.i;
#endif
    // increase flop counter
    fpar[10]=fpar[10]+2**n;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // left preconditioning
  if (lp) {
     // set status flag for a request for y = ML^{-1} r
     ipar[0]=3;
     // where does the source r start: w(:,1)
     ipar[7]=1;
     // where do we ask for the return value y: w(:,3)
     ipar[8]=2**n+1;
     // return flag, indicate that we want to re-enter at the fifth label 50
     ipar[9]=5;
     return;
  } // end if lp
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // we have computed y=ML^{-1}r

  // QMR update parameters
  label50:
    theta_old=theta;
    // write (6, '(A)') 'label 50'

    aalpha=fabs(tau);
    // |alpha|=0 or +/-oo ?
    if (BRKDN(&aalpha,ipar)) goto label900;
    i=1;
    j=1;
    if (lp) 
       theta=NRM(n,&w[0+2**n],&i)/tau;
    else
       theta=NRM(n,w,&i)/tau;

    // Givens rotation
    c2=1.0/(1+theta*theta);
    c=sqrt(c2);

    // norm of the precnd. quasi residual
    tau=tau*theta*c;
    fpar[4]=tau;

    // d=w(:,6)
    // p=w(:,2)

    // d = c^2 theta_old^2 d + c^2 alpha p
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    beta  =c2*theta_old*theta_old;
#else
    beta.r=c2*theta_old*theta_old;
    beta.i=0.0;
#endif
    SCAL(n,&beta,&w[0+5**n],&i);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    beta  =c2*alpha;
#else
    beta.r=c2*alpha.r;
    beta.i=c2*alpha.i;
#endif
    AXPY(n,&beta,&w[0+*n],&i,&w[0+5**n],&j);

    // x = x + d
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    alpha=1.0;
#else
    alpha.r=1.0; alpha.i=0.0;
#endif
    AXPY(n,&alpha,&w[0+5**n],&i,&w[0+4**n],&j);
    // increase flop counter
    fpar[10]=fpar[10]+7**n;
  // ------------------------------------------------------------------

     
  // ------------------------------------------------------------------
  // external stopping criterion
  if (ipar[2]==999) {
     // set status flag for a request for external stopping criterion
     ipar[0]=10;
     // x is stored in w(:,5)
     ipar[7]=4**n+1;
     // return result should be located in y=w(:,3)
     ipar[8]=2**n+1;
     // return flag, indicate that we want to re-enter at the sixth label 60
     ipar[9]=6;
     return;
  } // end if
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // return from external stopping criterion
  label60:
    if (ipar[2]==999) {
       if (ipar[10]==1) goto label900;
    }
    else if (ipar[2]<=2) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
       alpha  =1.0;
#else
       alpha.r=1.0;
       alpha.i=0.0;
#endif
       i=1;
       aalpha=FABS(alpha);
       if (STOPBIS(n,ipar,&i,fpar,w,&w[0+5**n],&aalpha)) goto label900;
    }
    else if (ipar[2]==3) {
       // delay update of ||x||
       if (delaycount>=delay) {
	  delaycount=1;
	  i=1;
	  nrmx=NRM(n,&w[0+4**n],&i);
	  // increase flop counter
	  ipar[10]=ipar[10]+2**n;
       }
       else
	  delaycount=delaycount+1;
       fpar[3]=adjust*fpar[0]*(nrmb+fpar[11]*nrmx);
       if ((ipar[6]>=ipar[5] && ipar[5]>0) || fpar[4]<=fpar[3]) {
	  // at first glance it looks as if we have succeeded
	  // to be save we check that the update has a small
	  // relative norm compared with x
	  //||d|| <= reltol ||x||/2
	  i=1;
	  nrmd=NRM(n,&w[0+5**n],&i);
	  // increase flop counter
	  ipar[10]=ipar[10]+2**n;
	  if ((ipar[6]>=ipar[5] && ipar[5]>0) || nrmd<=0.5*fpar[0]*nrmx) goto label900;
       } // end if
    } // end if
  // ------------------------------------------------------------------
 
  // ------------------------------------------------------------------
  // store old rho
  label310:
    rho_old=rho;
    fpar[7]=fpar[6];

    // request for right preconditioning t = MR^{-1} y

    // if right preconditioning is desired
    if (rp) { 
       // set status flag for a request for performing t = MR^{-1} y
       ipar[0]=5;
       if (lp) {
	  // y is stored in w(:,3)
	  ipar[7]=2**n+1;
	  // where do we ask for the return value t: w(:,4)
	  ipar[8]=3**n+1;
       }
       else {
	  // r is stored in w(:,1)
	  ipar[7]=1;
	  // where do we ask for the return value t: w(:,3)
	  ipar[8]=2**n+1;
       } // end if-else lp
       // return flag, indicate that we want to re-enter at the seventh label 70
       ipar[9]=7;
       return;
    } // end if rp
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // w(:,1) = r
  // w(:,2) = p
  // rho_new = t^* r
  label70:
    i=1;
    j=1;
    if ((lp!=0 && rp==0) || (lp==0 && rp!=0)) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
       rho=MYDISTDOT( n,&w[0+2**n],&i,w,&j);
#else
#ifdef _USE_MKL_
       MYDISTDOT(&rho,n,&w[0+2**n],&i,w,&j);
#else
       rho=MYDISTDOT( n,&w[0+2**n],&i,w,&j);
#endif 
#endif 
    }
    else {
       if (lp && rp) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	  rho=MYDISTDOT( n,&w[0+3**n],&i,w,&j);
#else
#ifdef _USE_MKL_
	  MYDISTDOT(&rho,n,&w[0+3**n],&i,w,&j);
#else
	  rho=MYDISTDOT( n,&w[0+3**n],&i,w,&j);
#endif 
#endif 
       }
       else {
	 // no preconditioning
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	 rho=MYDISTDOT( n,w,&i,w,&j);
#else
#ifdef _USE_MKL_
	 MYDISTDOT(&rho,n,w,&i,w,&j);
#else
	 rho=MYDISTDOT( n,w,&i,w,&j);
#endif 
#endif 
       }
    } // end if-else
    fpar[6]=FABS(rho);
    // increase flop counter
    fpar[10]=fpar[10]+2**n;
    // check for break down
    // |rho|=0 or +/-oo ?
    if (BRKDN(&fpar[6],ipar)) goto label900;

    // alpha = rho_new / rho_old
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
    alpha=rho/rho_old;
#else
    aalpha=1.0/FABS2(rho_old);
    beta.r= rho_old.r*aalpha;
    beta.i=-rho_old.i*aalpha;
    alpha.r=rho.r*beta.r-rho.i*beta.i;
    alpha.i=rho.r*beta.i+rho.i*beta.r;
#endif 

    // p = MR^{-1} ML^{-1} r  +  alpha p
    if ((lp!=0 && rp==0) || (lp==0 && rp!=0)) {
       for (i=0; i<*n; i++) {
	   // w(:,3) = ML^{-1} r   or   w(:,3) = MR^{-1} r
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	   w[i+*n]=w[i+2**n]+alpha*w[i+*n];
#else
	   beta.r=alpha.r*w[i+*n].r-alpha.i*w[i+*n].i;
	   beta.i=alpha.r*w[i+*n].i+alpha.i*w[i+*n].r;
	   w[i+*n].r=w[i+2**n].r+beta.r;
	   w[i+*n].i=w[i+2**n].i+beta.i;
#endif
       } // end for i
    }
    else {
      if (lp&&rp) {
	 for (i=0; i<*n; i++) {
	     // w(:,4) = MR^{-1} ML^{-1} r
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     w[i+*n]=w[i+3**n]+alpha*w[i+*n];
#else
	     beta.r=alpha.r*w[i+*n].r-alpha.i*w[i+*n].i;
	     beta.i=alpha.r*w[i+*n].i+alpha.i*w[i+*n].r;
	     w[i+*n].r=w[i+3**n].r+beta.r;
	     w[i+*n].i=w[i+3**n].i+beta.i;
#endif
	 } // end for i
      }
      else {
	 // no preconditioning
	 for (i=0; i<*n; i++) {
	     // w(:,1) = r
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
	     w[i+*n]=w[i]+alpha*w[i+*n];
#else
	     beta.r=alpha.r*w[i+*n].r-alpha.i*w[i+*n].i;
	     beta.i=alpha.r*w[i+*n].i+alpha.i*w[i+*n].r;
	     w[i+*n].r=w[i].r+beta.r;
	     w[i+*n].i=w[i].i+beta.i;
#endif
	 } // end for i
      } // end if-else lp&&rp
    } // end if-else lp != rp
    // increase flop counter
    fpar[10]=fpar[10]+2**n;


    // end of the iterations
    
    goto label300;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // some clean up job to do

  label900:
    if (ipar[11]!=0)
       ipar[0]=ipar[11];
    else if (ipar[0]>0) {
         if ((ipar[2]==999 && ipar[10]==1) || fpar[4]<=fpar[3]) 
	    ipar[0]=0;
         else if (ipar[6]>=ipar[5] && ipar[5]>0)
	    ipar[0]=-1;
         else
	    ipar[0]=-10;
    } // end if-else if
    if (fpar[2]>0.0 && fpar[5]>0.0 && ipar[6]>ipar[12])
       fpar[6]=log10(fpar[2]/fpar[5])/((double)(ipar[6]-ipar[12]));
    else
       fpar[6]=0.0;
    
    // failure
    if (ipar[0]<0) {
       i=1;
       j=1;
       COPY(n,&w[0+4**n],&i,sol,&j);
       return;
    } // end if
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // compute true norm of the residual

  // ------------------------------------------------------------------
  // request for matrix vector multiplication y = A x for finalization

  // set status flag for a request for performing y = A x
  i=ipar[0];
  ipar[0]=1;
  // where does the source x start: w(:,5)
  ipar[7]=4**n+1;
  // where do we ask for the return value y: w(:,3)
  ipar[8]=2**n+1;
  // return flag, indicate that we want to re-enter at eighth label 80
  ipar[9]=8;
  return;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // we have performed a matrix-vector multiplication
  // the result y=Ax is returned in w(:,3)
  // now compute the residual y=b-Ax

  // build residual y=b-Ax, store it in w(:,3)
  label80:
    ipar[0]=i;
    // write (6, '(A)') 'label 80'
    for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
        w[i+2**n]  =rhs[i]  -w[i+2**n];
#else
        w[i+2**n].r=rhs[i].r-w[i+2**n].r;
        w[i+2**n].i=rhs[i].i-w[i+2**n].i;
#endif
    } // end for i
    // flop counter
    fpar[10]=fpar[10]+*n;
      
    // ||b-Ax||
    i=1;
    aalpha=NRM(n,&w[0+2**n],&i);
    fpar[5]=aalpha;

    if (ipar[2]==3) {
       if (aalpha!=0.0) adjust=0.5*(fpar[4]/aalpha);
       // compute true backward error
       i=1;
       nrmx=NRM(n,&w[0+4**n],&i);
       fpar[3]=adjust*fpar[0]*(nrmb+fpar[11]*nrmx);
    } // end if

    // for the backward error the true residual has to be below
    // the threshold, 'adjust' is only used for safety
    if (ipar[2]==3 && fpar[5]<=fpar[0]*(nrmb+fpar[11]*nrmx)) goto label320;

    // if necessary continue QMR
    if ((ipar[6]<ipar[5] || ipar[5]<=0) && fpar[4]>fpar[3]) goto label310;

    label320:
      i=1;
      j=1;
      COPY(n,&w[0+4**n],&i,sol,&j);
      return;
} //end-of-SQMR
