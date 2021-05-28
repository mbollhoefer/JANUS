#include <ilupack_fortran.h>
#include <ilupackmacros.h>

      subroutine PCG(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      FLOAT  rhs(n), sol(n), w(n,*)
      REALS fpar(16)
c-----------------------------------------------------------------------
c     This is a implementation of the PRECONDITIONED Conjugate Gradient
c     (CG) method for solving linear system.
c
c-----------------------------------------------------------------------
c
c     n     size of the system
c
c     rhs   right hand side
c
c     sol   on entry, a given approximate solution
c           on exit, the new approximate solution
c
c     ipar  integer parameters
c   
c     fpar  real parameters
c 
c     w     work space
c
c-----------------------------------------------------------------------
c
c     ipar( 1):  status flag. depending on its values the algorithm is quit 
c                for different requests
c                0 initialization/return
c                1 matrix-vector multiplication y=Ax
c                2 transposed matrix-vector multiplication y=A^Tx
c                3 left preconditioning
c                4 left preconditioning transposed
c                5 right preconditioning
c                6 right preconditioning transposed
c               10 test external stopping routine
c               -1 too many iteration steps
c               -2 not enough work space
c               -3 algorithm breaks down
c                  for Arnoldi-type methods refer to ipar(12) for a
c                  more detailed description of the break-down
c
c     ipar( 2):  0 no preconditioning used
c                1 ONLY left preconditioning
c                2 ONLY right preconditioning
c                3 left AND right preconditioning
c
c     ipar( 3):  stopping criteria
c                  1  use residual norm as computed by cg
c                  2  use rel. tol. * ||b|| + abs. tol.
c                999  external stopping criterium, 
c
c     ipar( 4):  size of working array w
c
c     ipar( 5):  restart length (only used for GMRES,FOM,...
c
c     ipar( 6):  maximum number of iteration steps
c
c     ipar( 7):  counter for the number of steps (matrix-vector multiplications)
c
c     ipar( 8):  address of the source x when applying one of the operations
c                y in {Ax, A^Tx, Ml^{-1}x, Mr^{-1}x}
c
c     ipar( 9):  associated target address for y
c
c     ipar(10):  status flag. depending on its values the algorithm is 
c                reentered at different labels
c                0 initialization
c                1 return after INITIAL matrix-vector multiplication y=Ax
c                2 return after preconditioning
c                3 return after matrix-vector multiplication y=Ax

c     ipar(11):  flag for successful external stopping criterion
c                 0 failure
c                 1 success
c
c     ipar(12):  for Arnoldi-type methods: break down, reason
c		 -1  due to zero input vector
c		 -2  since input vector contains abnormal numbers
c		 -3  since input vector is a linear combination of others
c		 -4  since triangular system in GMRES/FOM/etc. has null rank
c                 
c     ipar(13):
c     ipar(14):
c     ipar(15):
c     ipar(16):
c
c
c     fpar( 1):  relative error tolerance
c
c     fpar( 2):  absolute error tolerance
c     fpar( 3):  initial bound of the residual
c     fpar( 4):  target bound for the residual
c     fpar( 5):  current error bound of the residual
c     fpar( 6):  current error bound of the delta x  (or the residual)
c     fpar( 7):  convergence rate, here: norm of the preconditioned residual squared
c     fpar( 8):
c     fpar( 9):
c     fpar(10):
c     fpar(11): flop counter
c     fpar(12):
c     fpar(13):
c     fpar(14):
c     fpar(15):
c     fpar(16):
c
c
c     fpar(7) is used here internally to store <r, r>.
c     w(:,1) -- residual vector
c     w(:,2) -- P, the conjugate direction
c     w(:,3) -- A P, matrix multiply the conjugate direction
c     w(:,4) -- temporary storage for results of preconditioning
c     w(:,5) -- change in the solution (sol) is stored here until
c               termination of this solver
c
c-----------------------------------------------------------------------
c     Authors:
c
c	Matthias Bollhoefer, TU Braunschweig
c
c     Date:
c
c       code modified from SPARSKIT-cg of Yousef Saad, adapted for PCG 
c       and the complex case
c       August 23, 2003. ILUPACK V1.0
c
c     Notice:
c
c	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
c	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
c
c     Availability:
c
c	This file is located at
c
c	http://ilupack.tu-bs.de/
c-----------------------------------------------------------------------
c     external functions used
c
      external DISTDOT,DISTDOT2,STOPBIS,BRKDN,BISINIT,NRM,COPY,AXPY
      integer iabs
      logical STOPBIS, BRKDN
      FLOAT  DISTDOT,DISTDOT2
      REALS NRM, SQRTFNAME
c
c     local variables
c
      integer i, d, nupos
      parameter(d=4)
      REALS alpha, rho, nrmx, nrmx0, nrmd, nrm0, 
     +      sumnu, nu(d) 
      FLOAT falpha
      logical lp,rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 30, 60, 70), ipar(10)
c
c     initialization
c
      call BISINIT(ipar,fpar,5*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
c     position inside array nu
      nupos=1
      nrmx0=-1.0

c     ------------------------------------------------------------------
c     request for matrix vector multiplication y=A*x for initialization
c
c     set status flag for a request for performing y=A*x
      ipar(1) = 1
c     where does the source x start: w(:,2)
      ipar(8) = n+1
c     where do we ask for the return value y: w(:,3)
      ipar(9) = ipar(8) + n
c     return flag, indicate that we want to re-enter at label 10
      ipar(10) = 1
c     copy the initial solution to the buffer for the source x
c     w(:,2)=sol
      i=1
      call COPY(n,sol,i,w(1,2),i)
      return
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     we have performed a matrix-vector multiplication
c     the result y=Ax is returned in w(:,3)
c     now compute the initial residual r=b-Ax
c     
c     increment step counter
 10   ipar(7) = ipar(7) + 1
      ipar(13) = 1
c     preparations if the relative error in the energy norm is required
      if (ipar(3).eq.3) then
         i=1
c        scalar product (b,x)
         falpha=DISTDOT(n,rhs,i,sol,i)
c        real part
         alpha=CTOD(falpha+CONJ(falpha))
c        2* real(b,x_0) - ||x_0||_A^2
         nrm0=alpha-CTOD(DISTDOT(n,sol,i,w(1,3),i))
         do i=1,d
            nu(i)=0.0
         end do
c        flop counter
         fpar(11) = fpar(11) + 4*n
      end if
c     build residual r=b-Ax, store it in w(:,1)
      do i = 1, n
         w(i,1) = rhs(i) - w(i,3)
      enddo
c     flop counter
      fpar(11) = fpar(11) + n
c     norm of the initial residual ||r||
      i=1
      fpar(3) = NRM(n,w,i)
c     flop counter
      fpar(11) = fpar(11) + 2*n
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     current norm, here: initial residual
      fpar(5) = fpar(3)
c     which stopping criterion did we select
c     case 1: reltol * ||b|| + abstol
      if (iabs(ipar(3)).eq.2) then
c        error bound fpar(1)=relative tolerance
c                    fpar(2)=absolute tolerance
c        model: reltol*||b||+abstol
         i=1
         fpar(4) = fpar(1)*SQRT(CTOD(DISTDOT2(n,rhs,i)))+fpar(2)
c        increase flop counter
         fpar(11) = fpar(11) + 2 * n
c     case 2: reltol * ||r|| + abstol
      else if (ipar(3).ne.999) then
c        error bound fpar(1)=relative tolerance
c                    fpar(2)=absolute tolerance
c        model: reltol*||r||+abstol, where r is the residual
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
c        fpar(4) not used if ipar(3)=3
      endif
c     ------------------------------------------------------------------



c     ------------------------------------------------------------------
c     ---------------------     START MAIN LOOP     --------------------
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     if preconditioning z=M^{-1}r is desired 
c
 20   if (lp .or. rp) then
c        set status flag for a request for preconditioning
         if (lp) then
            ipar(1) = 3
         else
            ipar(1) = 5
         endif
c        where does the source r start: w(:,1)
         ipar(8) = 1
c        target z for preconditioning z=M{-1}r: w(:,4)
         ipar(9) = 3*n+1
c        return flag, indicate to return to label 30 when re-entering 
c        the subroutine
         ipar(10) = 2
         return
      else
c        use w(:,4)=w(:,1)
         i=1
         call COPY(n,w,i,w(1,4),i)
      endif
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     compute scalar product rho=r^Hz of between the residual r and the 
c     preconditioned residual z. rho must be real and positive
 30   i=1
      rho=CTOD(DISTDOT(n,w(1,1),i,w(1,4),i))
c     increment flop counter
      fpar(11)= fpar(11) + 2*n
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     init/update conjugate search direction
      if (ipar(7).eq.1) then
c        initial conjugate direction stored in w(:,2)
c        p=z
         i=1
         call COPY(n,w(1,4),i,w(1,2),i)
      else
c        parameter for update of the search direction
c        r_new^Hz_new / r_old^Hz_old
         alpha = rho / fpar(7)
c        p = z + alpha p
         do i = 1, n
            w(i,2) = w(i,4) + alpha * w(i,2)
         enddo
c        increment flop counter
         fpar(11) = fpar(11) + 2*n
      endif
c     use fpar(7) to store the actual rho
      fpar(7) = rho
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     before iteration can continue, we need to compute A * p, where
c     p=w(:,2) is the conjugate search direction
c     set status flag as request for matrix-vector multiplication w=Ap
 40   ipar(1) = 1
c     source address: w(:,2), conjugate search direction p
      ipar(8) = n + 1
c     address of the target w: w(:,3)
      ipar(9) = 2*n+1
c     return flag, indicate to return to label 60 when re-entering the
c     subroutine
      ipar(10) = 3
      return
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     continueing with the iterations. 
c
c     advance counter for the matrix-vector multiplications
 60   ipar(7) = ipar(7) + 1
c     compute alpha=p^HAp
      i=1
      alpha = CTOD(DISTDOT(n,w(1,2),i,w(1,3),i))
c     increment flop counter
      fpar(11) = fpar(11) + 2*n
c     alpha=0 or +/-oo ?
      if (BRKDN(alpha,ipar)) goto 900
c     parameter to update approx. solution and residual
c     dx = dx + (r^Hz/p^HAp)  p, dx: sol+w(:,5) approx. sol.
c     r  = r  - (r^Hz/p^HAp) Ap,  r: w(:,1)     residual
c                                 p: w(:,2)     conjg. search direction
c                                Ap: w(:,3)     A * conjg. search direction
c     dx = dx + alpha  p
      falpha = rho/alpha
      if (ipar(3).eq.3) then
c        additional contribution to the energy norm estimate
         nu(nupos)=rho*rho/alpha
c        compute nu_j,d
         sumnu=0.0
         do i=1,d
            sumnu=sumnu+nu(i)
         end do
         fpar(5)=SQRT(abs(sumnu))
c        update nu_0,j + 2* real(b,x)-||x_0||_A^2
         nrm0=nrm0+nu(nupos)
         nupos=nupos+1
         if (nupos.gt.d) nupos=1
      end if
      i=1
      call AXPY(n,falpha,w(1,2),i,w(1,5),i)
c     r = r   - alpha Ap
      falpha=-falpha
      call AXPY(n,falpha,w(1,3),i,w(1,1),i)
      alpha=CTOD(-falpha)
c     increment flop counter
      fpar(11) = fpar(11) + 4*n
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     are we ready to terminate ?
c
c     external stopping criterion
      if (ipar(3).eq.999) then
c        set status flag for a request for external stopping criterion
         ipar(1) = 10
c        first address w(:,5) (approx. solution dx)
         ipar(8) = 4*n + 1
c        second address w(:,4), some target
         ipar(9) = 3*n + 1
c        return flag, indicate to return to label 70 when re-entering the
c        subroutine
         ipar(10) = 4
         return
      endif
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     check stopping criteria
 70   if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
c        w(:,1) residual
c        w(:,2) conjg. search direction
c        alpha*w(:,2) ~ delta x change in the solution
      else if (ipar(3).eq.3 .and. ipar(7).gt.d) then 
         if ((ipar(7).gt.d .and. sumnu.lt.fpar(1)*fpar(1)*(nrm0+sumnu))
     +       .or. ipar(7).ge.ipar(6) .or. alpha.eq.0.0) then
            fpar(4)=fpar(1)*SQRT(nrm0+sumnu)
c           at first glance it looks as if we have succeeded
c           but to be save we check that the update has a small
c           relative norm compared with x
c           ||p|| <= reltol ||x||/2
c           for ||x|| we use ||x0||+||dx|| as estimate
            nrmd=NRM(n,w(1,2),1)
            nrmx=NRM(n,w(1,5),1)
            if (nrmx0.lt.0.0) nrmx0=NRM(n,sol,1)
c           increase flop counter
            ipar(11) = ipar(11) + 4*n
            if ((ipar(7).ge.ipar(6)) .or.
     +          nrmd.le.0.5*fpar(1)*(nrmx+nrmx0)) goto 900
         else
         end if
      else if (ipar(3).ne.3 .or. (ipar(3).eq.3 .and. ipar(7).le.d)) then
         if (STOPBIS(n,ipar,1,fpar,w,w(1,2),alpha)) goto 900
      endif
c     ------------------------------------------------------------------


c     continue the iterations
      goto 20
c     ------------------------------------------------------------------
c     ----------------------     END MAIN LOOP     ---------------------
c     ------------------------------------------------------------------



c     ------------------------------------------------------------------
c     clean up 
c
c     sol = sol + w(:,5)
 900  call TIDYCG(n,ipar,fpar,sol,w(1,5))
c
      return
      end
c-----end-of-pcg
