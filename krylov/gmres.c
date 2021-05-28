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




void GMRES(integer *n, FLOAT *rhs, FLOAT *sol, integer *ipar, REALS *fpar, FLOAT *w)
{
  /* $Id: gmres.c 6094 2020-03-28 16:17:50Z bolle $

    -----------------------------------------------------------------------
         This a version of GMRES implemented with reverse communication.
         It is a simple restart version of the GMRES algorithm.
    
         ipar(1) control flag for external operation
                  1  y=Ax       matrix-vector multiplication
                  2  y=A^Tx     transposed matrix-vector multiplication
                  3  z=ML^{-1}r left preconditioning
                  4  z=ML^{-T}r transposed left preconditioning
                  5  z=MR^{-1}r right preconditioning
                  6  z=MR^{-T}r transposed right preconditioning
                 10  external stopping criterion
    
         ipar(3) stopping criterion, model: r_k=b-A*x_k
                  1 ||ML^{-1}*r_k||  <=  reltol*||r_0||+abstol
                  2 ||ML^{-1}*r_k||  <=  reltol*||b||  +abstol
                  3 backward error
                    ||r_k||          <=  reltol*(||b||+||A||*||x_k||)
    
         ipar(5) == the dimension of the Krylov subspace
         after every ipar(5) iterations, the GMRES will restart with
         the updated solution and recomputed residual vector.
    
         the space of the `w' is used as follows:
         (1) the basis for the Krylov subspace, size n*(m+1);
         (2) the Hessenberg matrix, only the upper triangular
         portion of the matrix is stored, size (m+1)*m/2 + 1
         (3) three vectors, all are of size m, they are
         the cosine and sine of the Givens rotations, the third one holds
         the residuals, it is of size m+1.
    
         TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
         Note: m == ipar(5). The default value for this is 15 if
         ipar(5) <= 1.
    
         code taken from SPARSKIT of Yousef Saad.
         adapted by Matthias Bollhoefer for the complex case
    -----------------------------------------------------------------------
  */

  static REALS one=1.0;

  /*  
c     local variables, ptr and p2 are temporary pointers,
c     hess points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
  */
  static integer i,j,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn;
  static REALS alpha, nrmb,nrmx,alpha0,adjust,deltax,c;
  static FLOAT  beta, s;
  static integer lp, rp, flag;

  //  check the status of the call

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

  
  // initialization

  // restart length for GMRES(m)
  if (ipar[4]<=1) 
     m=15;
  else
     m=ipar[4];

  // size of the matrix holding the m+1 Arnoldi vectors
  idx=*n*(m+1);
  // space for n additional elements
  hess=idx+*n;
  // space for the upper Hessenberg matrix
  vc=hess+(m+1)*m/2+1;
  // space for m cosine elements
  vs=vc+m;
  // space for m sine elements
  vrn=vs+m;
  // space for m+1 elements of the right hand side from ||alpha e_1 -Hy||
  i=vrn+m+1;
  j=1;
  BISINIT(ipar,fpar,&i,&j,&lp,&rp,w);
  if (ipar[0]<0) return;
  // adjust residual and preconditioned residual
  adjust=1.0;

  // ------------------------------------------------------------------
  // request for matrix vector multiplication y=A*x for initialization

  // set status flag for a request for performing y=A*x
  label100:
    ipar[0]=1;
    // where does the source x start: w(:,2)
    ipar[7]=*n+1;
    // where do we ask for the return value y: w(:,1)
    ipar[8]=1;
    // return flag, indicate that we want to re-enter at the first label 10
    ipar[9]=1;
    // reset counter for a restarted GMRES sequence      
    k=0;

    // copy the initial solution to the buffer for the source x
    // w(:,2)=sol
    j=1;
    COPY(n,sol,&j,&w[*n],&j);
    return;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // return from matrix-vector multiplication
  // increment step counter
  label10:
    ipar[ 6]=ipar[ 6]+1;
    ipar[12]=ipar[12]+1;
    if (lp) {
       // build residual r=b-Ax, store it in w(:,2), this is the place
       // for the source for left preconditioning
       for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   w[*n+i]=rhs[i]-w[i];
#else
	   w[*n+i].r=rhs[i].r-w[i].r;
	   w[*n+i].i=rhs[i].i-w[i].i;
#endif
       } // end for i
       // alpha0 = ||r||
       i=1;
       alpha0=NRM(n,&w[*n],&i);

       // flop counter
       fpar[10]=fpar[10]+3**n;

       // set status flag for a request for z=ML^{-1}*r
       ipar[0]=3;
       // return flag, indicate that we want to re-enter at the second label 20
       ipar[9]=2;
       return;
    } // end if
    else { // !lp
       // build residual r=b-Ax, rewrite it to w(:,1), this is the place where
       // left preconditioning would have returned the result
       for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   w[i]=rhs[i]-w[i];
#else
	   w[i].r=rhs[i].r-w[i].r;
	   w[i].i=rhs[i].i-w[i].i;
#endif
       } // end for i
       // alpha0 = ||r||
       i=1;
       alpha0=NRM(n,w,&i);

       // flop counter
       fpar[10]=fpar[10]+3**n;
    } // end if-else lp
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // alpha = ||ML^{-1} r||
  label20:
    i=1;
    if (lp) {
       alpha=NRM(n,w,&i);
       // flop counter
       fpar[10]=fpar[10]+2**n;
    } // end if
    else // !lp
       alpha=alpha0;
    // check step counter and stopping criterion
    if (ipar[6]==1 && ipar[2]!=999) {
       // which stopping criterion did we select
       // case 1: reltol * ||b|| + abstol
       if (IABS(ipar[2])==2) {
	  // error bound fpar(1)=relative tolerance
	  //             fpar(2)=absolute tolerance
	  // model: reltol*||b||+abstol
	  nrmb=NRM(n,rhs,&i);
	  fpar[3]=fpar[0]*nrmb+fpar[1];
	  // flop counter
	  fpar[10]=fpar[10]+2**n;
       }
       else if (IABS(ipar[2])==1) {
	  // case 2: reltol * ||ML^{-1}r_0|| + abstol
	  // error bound fpar(1)=relative tolerance
	  //             fpar(2)=absolute tolerance
	  // model: reltol*||ML^{-1}r_0||+abstol, where r is the residual
	  fpar[3]=fpar[0]*alpha+fpar[1];
       }
       else {
	  // case 3: ||r_k|| <= reltol * (||b|| + ||A||*||x_k||)
	  nrmb=NRM(n,rhs,&i);
	  nrmx=NRM(n,sol,&i);
	  deltax=nrmx;
	  // fpar[3]=fpar[0]*(nrmb+fpar[11]*nrmx)
	  fpar[3]=fpar[0]*nrmb;
	  // flop counter
	  fpar[10]=fpar[10]+4**n;
       } // end if-else if-else
       // ||ML^{-1}r_0 ||, initial residual
       fpar[2]=alpha;
    }
    else {
       if (IABS(ipar[2])==3) {
	  if (deltax<=fpar[0]*nrmx)
	     fpar[3]=fpar[0]*(nrmb+fpar[11]*nrmx);
	  else
	    fpar[3]=fpar[0]*nrmb;
       } // end if
    } // end if-else
    // ||ML^{-1}r ||, actual preconditioned residual
    fpar[4]=alpha;
    // leading component of the right hand side in 
    // ||ML^{-1}(b-Ax)||=||alpha e_1 - Hy||
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    w[vrn]=alpha;
#else
    w[vrn].r=alpha; w[vrn].i=0.0;
#endif

    // ||ML^{-1}r||<= target     &   residual norm 
    if (alpha<=fpar[3] && ipar[2]>=0 && ipar[2]!=999 && ipar[2]!=3) {
       // exit
       ipar[0]=0;
       // norm ||ML^{-1}(b-Ax)|| on exit
       fpar[5]=alpha;
       goto label300;
    } // end if 
    else if (alpha0<=fpar[3] && ipar[2]==3) {
       // exit
       ipar[0]=0;
       // norm ||(b-Ax)|| on exit
       fpar[5]=alpha0;
       goto label300;
    } // end if-else if
    // reciprocal for V(:,1)=ML^{-1}r/alpha
    alpha=one/alpha;

    // normalize w
    j=1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    beta=alpha;
#else
    beta.r=alpha; beta.i=0.0;
#endif
    SCAL(n,&beta,w,&j);

    // flop counter
    fpar[10]=fpar[10]+*n;


  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  //                         MAIN LOOP
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  // request for (1) right preconditioning
  //             (2) matrix vector multiplication
  //             (3) left preconditioning

  // increase number of steps until restart
  label110:
    k=k+1;

    // w(:,k+1)=ML^{-1}AMR^{-1}w(:,k)

    // right preconditioning
    if (rp) {
       // request for right preconditioning v=MR^{-1}u
       ipar[0]=5;
       // source u=w(:,k)
       ipar[7]=k**n-*n+1;
       if (lp)
	  // result v=w(:,k+1)
	  ipar[8]=k**n+1;
       else
	  // result z=w(:,m+2), m+1 vectors are needed for restart
	  ipar[8]=idx+1;
       // return flag, indicate that we want to re-enter at the third label 30
       ipar[9]=3;
       return;
    } // end if
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // request for x=Av
  label30:
    ipar[0]=1;
    // write (6,'(A)') 'return from prec., apply mat-vec'
    if (rp) 
       // if right preconditioning has been applied, then source v is located
       // at the position of the result v from v=MR^{-1}u
       ipar[7]=ipar[8];
    else
       // source v=w(:,k)
       ipar[7]=(k-1)**n+1;
    // if left preconditioning is also desired 
    if (lp)
       // result x=w(:,m+2), m+1 vectors are needed for restart
       ipar[8]=idx+1;
    else
       // no left prec., result x=w(:,k+1), space for the next Arnoldi vector
       ipar[8]=1+k**n;
    // return flag, indicate that we want to re-enter at the fourth label 40
    ipar[9]=4;
    return;
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // left preconditioning y=ML^{-1}x
  label40:
    if (lp) {
       // request for left preconditioning y=ML^{-1}x
       ipar[0]=3;
       // source x is obtained from previous step x=Av
       ipar[7]=ipar[8];
       // result y is returned in y=w(:,k+1)
       ipar[8]=k**n+1;
       // return flag, indicate that we want to re-enter at the fifth label 50
       ipar[9]=5;
       return;
    } // end if lp
    // write (6,'(A)') 'return from mat-vec., compute Gram-Schmidt'
  // ------------------------------------------------------------------

  // Modified Gram-Schmidt orthogonalization procedure
  // temporary pointer 'ptr' is pointing to the current column of the
  // Hessenberg matrix. 'p2' points to the new basis vector

  // increment step counter
  label50:
    ipar[6]=ipar[6]+1;
    // space to column of H(:,k+1), the previous columns are already transformed
    // to triangular form using Givens rotations
    ptr=k*(k-1)/2+hess;
    p2=ipar[8];
    // reorthogonalize w(:,k+1) against its predecessors using MGS with
    // selective reorthogonalization
    flag=0;
    k=k+1;
    MGSRO(&flag,n,n,&k,&k,&fpar[10],w,&w[ptr],&ipar[11]);
    k=k-1;
    // write (6,'(A,I4)') 'Gram-Schmidt completed, ipar(12)=',ipar(12)
    if (ipar[11]<0) goto label200;

    // apply previous Givens rotations and generate a new one to eliminate
    // the subdiagonal element.

    // position of H(:,k)
    p2=ptr+1;
    // H(1:k,k) is transformed using of old rotations
    ii=1;
    for (i=0; i<k-1; i++) {
        // H(i:i+1,k) is rotated
        ptr=p2;
        p2=p2+1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	alpha=w[vc+i];
#else
	alpha=w[vc+i].r;
#endif
	ROT(&ii,&w[ptr-1],&ii,&w[p2-1],&ii,&alpha,&w[vs+i]);
    } // end for i
    // write (6,'(A)') 'last column of H updated'
    // eliminate H(k:k+1,k) by Givens rotation
    ROTG(&w[p2-1],&w[p2],&c,&s);
    // store c and s
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    w[vc+k-1]=c;
#else
    w[vc+k-1].r=c; w[vc+k-1].i=0.0;
#endif
    w[vs+k-1]=s;
    // update right hand side
    // [alpha; 0] -> [c; -s']*alpha
    p2=vrn+k;
    
    // beta=-CONJ(s)*w[p2-1];
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    beta=-s*w[p2-1];
    w[p2-1]=c*w[p2-1];
#else
    beta.r=-(s.r*w[p2-1].r+s.i*w[p2-1].i);
    beta.i=-(s.r*w[p2-1].i-s.i*w[p2-1].r);
    w[p2-1].r=c*w[p2-1].r;
    w[p2-1].i=c*w[p2-1].i;
#endif

    w[p2]  =beta;
    // write (6,'(A)') 'end Arnoldi'

    // end of one Arnoldi iteration, alpha will store the estimated
    // residual norm at current stage

    // flop counter
    fpar[10]=fpar[10]+6*k+2;
    // new residual norm
    alpha=FABS(beta);
    fpar[4]=alpha;

    // k<m:                          within a restarted sequence
    // ipar(3)>=0 & alpha<=fpar(4):  ||ML^{-1}r||<=threshold (relative and/or absolute)
    // ipar(6)<=0:                   no max. its.
    // ipar(7)<ipar(6):              number of iterations less than max. its.
    if (k<m && !(ipar[2]>=0 && alpha<=adjust*fpar[3]) && (ipar[5]<=0 || ipar[6]<ipar[5])) goto label110;
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  //                         END MAIN LOOP
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------


  // update the approximate solution, first solve the upper triangular
  // system, temporary pointer ptr points to the Hessenberg matrix,
  // p2 points to the right-hand-side (also the solution) of the system.

  // ------------------------------------------------------------------
  label200:
    ptr=hess+k*(k+1)/2;
    p2=vrn+k;
    if (FABS(w[ptr-1])==0.0) {
      // if the diagonal elements of the last column is zero, reduce k by 1
      // so that a smaller triangular system is solved [It should only
      // happen when the matrix is singular, and at most once!]

      k=k-1;
      if (k>0) 
	 goto label200;
      else {
	 ipar[0]=-3;
	 ipar[11]=-4;
	 goto label300;
      } // end if-else
    } // end if
    // H -> Givens -> R
    // backward solve R*y=w(vrn+(1:k))
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
    w[p2-1]=w[p2-1]/w[ptr-1];
#else
    c=w[ptr-1].r*w[ptr-1].r+w[ptr-1].i*w[ptr-1].i;
    s.r= w[ptr-1].r/c;
    s.i=-w[ptr-1].i/c;
    beta.r=w[p2-1].r*s.r-w[p2-1].i*s.i;
    beta.i=w[p2-1].r*s.i+w[p2-1].i*s.r;
    w[p2-1]=beta;
#endif
    j=1;
    for (i=k-1; i>=1; i--) {
        ptr=ptr-i-1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	w[p2-1]=-w[p2-1];
#else
	w[p2-1].r=-w[p2-1].r;
	w[p2-1].i=-w[p2-1].i;
#endif
	AXPY(&i,&w[p2-1],&w[ptr],&j,&w[vrn],&j);
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	w[p2-1]=-w[p2-1];
#else
	w[p2-1].r=-w[p2-1].r;
	w[p2-1].i=-w[p2-1].i;
#endif

	p2=p2-1;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	w[p2-1]=w[p2-1]/w[ptr-1];
#else
	c=w[ptr-1].r*w[ptr-1].r+w[ptr-1].i*w[ptr-1].i;
	s.r= w[ptr-1].r/c;
	s.i=-w[ptr-1].i/c;
	beta.r=w[p2-1].r*s.r-w[p2-1].i*s.i;
	beta.i=w[p2-1].r*s.i+w[p2-1].i*s.r;
	w[p2-1]=beta;
#endif
    } // end for i
  // ------------------------------------------------------------------

  // update approximate solution
     
  // ------------------------------------------------------------------
  // first compute w(:,1)=Vy
  SCAL(n,&w[p2-1],w,&j);

  for (i=1; i<=k-1; i++) {
      ptr=i**n;
      p2=p2+1;

      AXPY(n,&w[p2-1],&w[ptr],&j,w,&j);
  } // end for i
  // flop counter
  fpar[10]=fpar[10]+2*k**n-*n+k*(k+1);
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // replace w(:,1)=Vy by MR^{-1} w(:,1)
  if (rp) {
     // request for right preconditioning
     ipar[0]=5;
     // source w(:,1)
     ipar[7]=1;
     // result w(:,m+2)
     ipar[8]=idx+1;
     // return to sixth label 60
     ipar[9]=6;
     return;
  } // end if rp
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // update approximate solution
  label60:
    if (rp) {
       for (i=0; i<*n; i++) {
	   // update solution using MR^{-1} Vy
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   sol[i]=sol[i]+w[idx+i];
#else
	   sol[i].r=sol[i].r+w[idx+i].r;
	   sol[i].i=sol[i].i+w[idx+i].i;
#endif
       } // end for i
       if (IABS(ipar[2])==3) {
	  i=1;
	  deltax=NRM(n,&w[idx-1],&i);
       } // end if
    } // end if 
    else { // !rp
       for (i=0; i<*n; i++) {
	   // without right preconditioning we can directly update sol using w(:,1)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   sol[i]=sol[i]+w[i];
#else
	   sol[i].r=sol[i].r+w[i].r;
	   sol[i].i=sol[i].i+w[i].i;
#endif
       } // end for i
       if (IABS(ipar[2])==3) {
	  i=1;
	  deltax=NRM(n,w,&i);
       } // end if
    } // end if-else rp
    // increase flop counter
    fpar[10]=fpar[10]+*n;
  // ------------------------------------------------------------------

  // process the complete stopping criteria
  if (ipar[2]==999) {
     ipar[0]=10;
     ipar[7]=-1;
     ipar[8]=idx+1;
     ipar[9]=7;
     return;
  } // end if
  else if (ipar[2]<0) {
    if (ipar[6]<=m+1) {
       fpar[2]=FABS(w[vrn]);
       if (ipar[2]==-1) fpar[3]=fpar[0]*fpar[2]+fpar[1];
    } // end if
    fpar[5]=FABS(w[vrn+k-1]);
  }
  else if (ipar[2]==1 || ipar[2]==2) {
     // store norm of the preconditioned residual
     fpar[5]=fpar[4];
  } 
  else {
     i=1;
     // update right hand side for the backward error
     nrmx=NRM(n,sol,&i);
     if (deltax<=fpar[0]*nrmx) 
        fpar[3]=fpar[0]*(nrmb+fpar[11]*nrmx);
     else
        fpar[3]=fpar[0]*nrmb;
     // fpar[3]=fpar[0]*(nrmb+fpar[11]*nrmx)
     
     // ------------------------------------------------------------------
     // compute exact residual b-Ax

     // set status flag for a request for performing y=A*x
     ipar[0]=1;
     // where does the source x start: w(:,2)
     ipar[7]=*n+1;
     // where do we ask for the return value y: w(:,1)
     ipar[8]=1;
     // return flag, indicate that we want to re-enter at the eighth label 80
     ipar[9]=8;
     // reset counter for a restarted GMRES sequence      
     k=0;

     // copy the initial solution to the buffer for the source x
     // w(:,2)=sol
     j=1;
     COPY(n,sol,&j,&w[*n],&j);
     return;
     // ------------------------------------------------------------------
  } // end if-else if-else


  // do we need to restart ?

  label70:
    if (ipar[11]!=0) {
       ipar[0]=-3;
       goto label300;
    } // end if
    if ((ipar[6]<ipar[5] || ipar[5]<=0) &&
	((ipar[2]==999 && ipar[10]==0) || (ipar[2]!=999 && fpar[5]>fpar[3]))) goto label100;

  // termination, set error code, compute convergence rate

  label310:
    if (ipar[0]>0) {
       if (ipar[2]==999 && ipar[10]==1)
	  ipar[0]=0;
       else if (ipar[2]!=999 && fpar[5]<=fpar[3])
	  ipar[0]=0;
       else if (ipar[6]>=ipar[5] && ipar[5]>0)
	  ipar[0]=-1;
       else
	  ipar[0]=-10;
    } // end if
       
  label300:
    if (fpar[2]!=0.0 && fpar[5]!=0.0 && ipar[6]>ipar[12]) 
       fpar[6]=log10(fpar[2]/fpar[5])/((REALS)(ipar[6]-ipar[12]));
    else
       fpar[6]=0.0;
    return;
  // ------------------------------------------------------------------


  // ------------------------------------------------------------------
  // return from matrix-vector multiplication
  label80:
    if (lp) {
       // build residual r=b-Ax, store it in w(:,2), this is the place
       // for the source for left preconditioning
       for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   w[*n+i]=rhs[i]-w[i];
#else
	   w[*n+i].r=rhs[i].r-w[i].r;
	   w[*n+i].i=rhs[i].i-w[i].i;
#endif
       } // end for i
       // alpha0 = ||r||
       i=1;
       alpha0=NRM(n,&w[*n],&i);

       // flop counter
       fpar[10]=fpar[10]+3**n;
       
       // set status flag for a request for z=ML^{-1}*r
       ipar[0]=3;
       // return flag, indicate that we want to re-enter at the second label 20
       ipar[9]=2;
    } // end if
    else {
       // build residual r=b-Ax, rewrite it to w(:,1), this is the place where
       // left preconditioning would have returned the result
       for (i=0; i<*n; i++) {
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_		     
	   w[i]=rhs[i]-w[i];
#else
	   w[i].r=rhs[i].r-w[i].r;
	   w[i].i=rhs[i].i-w[i].i;
#endif
       } // end for i
       // alpha0 = ||r||
       i=1;
       alpha0=NRM(n,w,&i);
       
       // flop counter
       fpar[10]=fpar[10]+3**n;
    } // end if-else
  // ------------------------------------------------------------------

  if (ipar[11]!=0) {
     ipar[0]=-3;
     goto label300;
  } // end if
  // store norm of the exact residual
  fpar[5]=alpha0;
  // adjust relation between the original and the preconditioned residual
  // only the preconditioned residual can be easliy measured during a single
  // restart sequence. For safety, add a factor 0.5
  adjust=0.5*(fpar[4]/alpha0);
  if ((ipar[6]<ipar[5] || ipar[5]<=0) && fpar[5]>fpar[3]) {
     // restart GMRES

     // increment step counter
     ipar[ 6]=ipar[ 6]+1;
     ipar[12]=ipar[12]+1;
     if (lp) 
        // call for left preconditioning to re-enter restarted GMRES
        return;
     else
        // re-enter restarted GMRES without left preconditioning
        goto label20;
  } // end if
  goto label310;
} // end-of-gmres
