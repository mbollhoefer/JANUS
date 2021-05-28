#include <ilupack_fortran.h>
#include <ilupackmacros.h>



      subroutine GMRES(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      FLOAT  rhs(n), sol(n), w(*)
      REALS   fpar(16)
c-----------------------------------------------------------------------
c     This a version of GMRES implemented with reverse communication.
c     It is a simple restart version of the GMRES algorithm.
c
c     ipar(1) control flag for external operation
c              1  y=Ax       matrix-vector multiplication
c              2  y=A^Tx     transposed matrix-vector multiplication
c              3  z=ML^{-1}r left preconditioning
c              4  z=ML^{-T}r transposed left preconditioning
c              5  z=MR^{-1}r right preconditioning
c              6  z=MR^{-T}r transposed right preconditioning
c             10  external stopping criterion
c
c     ipar(3) stopping criterion, model: r_k=b-A*x_k
c              1 ||ML^{-1}*r_k||  <=  reltol*||r_0||+abstol
c              2 ||ML^{-1}*r_k||  <=  reltol*||b||  +abstol
c              3 backward error
c                ||r_k||          <=  reltol*(||b||+||A||*||x_k||)
c
c     ipar(5) == the dimension of the Krylov subspace
c     after every ipar(5) iterations, the GMRES will restart with
c     the updated solution and recomputed residual vector.
c
c     the space of the `w' is used as follows:
c     (1) the basis for the Krylov subspace, size n*(m+1);
c     (2) the Hessenberg matrix, only the upper triangular
c     portion of the matrix is stored, size (m+1)*m/2 + 1
c     (3) three vectors, all are of size m, they are
c     the cosine and sine of the Givens rotations, the third one holds
c     the residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c
c     code taken from SPARSKIT of Yousef Saad.
c     adapted by Matthias Bollhoefer for the complex case
c-----------------------------------------------------------------------
c     external functions used
c
      external DISTDOT, DISTDOT2, AXPY, COPY, NRM, BISINIT, 
     +         ROT,ROTG, SCAL
      FLOAT DISTDOT, DISTDOT2
      integer iabs
c
      REALS ABSFNAME, LOGTENFNAME, one, zero, NRM
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_
      parameter(one=1.0e0, zero=0.0e0)
#else
      parameter(one=1.0D0, zero=0.0D0)
#endif
c
c     local variables, ptr and p2 are temporary pointers,
c     hess points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,j,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      REALS alpha, nrmb,nrmx,alpha0,adjust,deltax
      FLOAT  beta, c, s
      logical lp, rp
      save

c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70, 80) ipar(10)
c
c     initialization
c
c     restart length for GMRES(m)
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
c     size of the matrix holding the m+1 Arnoldi vectors
      idx = n * (m+1)
c     space for n additional elements
      hess = idx + n
c     space for the upper Hessenberg matrix
      vc = hess + (m+1) * m / 2 + 1
c     space for m cosine elements
      vs = vc + m
c     space for m sine elements
      vrn = vs + m
c     space for m+1 elements of the right hand side from ||alpha e_1 -Hy||
      i = vrn + m + 1
      call BISINIT(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
c     adjust residual and preconditioned residual
      adjust=1.0

c     ------------------------------------------------------------------
c     request for matrix vector multiplication y=A*x for initialization
c
c     set status flag for a request for performing y=A*x
 100  ipar(1) = 1
c     where does the source x start: w(:,2)
      ipar(8) = n+1
c     where do we ask for the return value y: w(:,1)
      ipar(9) = 1
c     return flag, indicate that we want to re-enter at the first label 10
      ipar(10) = 1
c     reset counter for a restarted GMRES sequence      
      k = 0

c     copy the initial solution to the buffer for the source x
c     w(:,2)=sol
      j=1
      call COPY(n,sol,j,w(n+1),j)
      return
c     ------------------------------------------------------------------

c     ------------------------------------------------------------------
c     return from matrix-vector multiplication
c     increment step counter
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
c        build residual r=b-Ax, store it in w(:,2), this is the place
c        for the source for left preconditioning
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         end do
c        alpha0 = ||r||
         i=1
         alpha0 = NRM(n,w(n+1),i)

c        flop counter
         fpar(11) = fpar(11) + 3*n

c        set status flag for a request for z=ML^{-1}*r
         ipar(1) = 3
c        return flag, indicate that we want to re-enter at the second label 20
         ipar(10) = 2
         return
      else
c        build residual r=b-Ax, rewrite it to w(:,1), this is the place where
c        left preconditioning would have returned the result
         do i = 1, n
            w(i) = rhs(i) - w(i)
         end do
c        alpha0 = ||r||
         i=1
         alpha0 = NRM(n,w,i)

c        flop counter
         fpar(11) = fpar(11) + 3*n
      endif
c     ------------------------------------------------------------------

c     ------------------------------------------------------------------
c     alpha = ||ML^{-1} r||
 20   i=1
      if (lp) then
         alpha = NRM(n,w,i)
c        flop counter
         fpar(11) = fpar(11) + 2*n
      else
         alpha=alpha0
      end if
c     check step counter and stopping criterion
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
c        which stopping criterion did we select
c        case 1: reltol * ||b|| + abstol
         if (iabs(ipar(3)).eq.2) then
c           error bound fpar(1)=relative tolerance
c                       fpar(2)=absolute tolerance
c           model: reltol*||b||+abstol
            nrmb=NRM(n,rhs,i)
            fpar(4) = fpar(1)*nrmb+fpar(2)
c           flop counter
            fpar(11) = fpar(11) + 2*n
         elseif (iabs(ipar(3)).eq.1) then
c        case 2: reltol * ||ML^{-1}r_0|| + abstol
c           error bound fpar(1)=relative tolerance
c                       fpar(2)=absolute tolerance
c           model: reltol*||ML^{-1}r_0||+abstol, where r is the residual
            fpar(4) = fpar(1) * alpha + fpar(2)
         else
c        case 3: ||r_k|| <= reltol * (||b|| + ||A||*||x_k||)
            nrmb=NRM(n,rhs,i)
            nrmx=NRM(n,sol,i)
            deltax=nrmx
c            fpar(4)=fpar(1)*(nrmb+fpar(12)*nrmx)
            fpar(4)=fpar(1)*nrmb
c           flop counter
            fpar(11) = fpar(11) + 4*n
         end if
c        ||ML^{-1}r_0 ||, initial residual
         fpar(3) = alpha
      else
         if (iabs(ipar(3)).eq.3) then
            if (deltax.le.fpar(1)*nrmx) then
               fpar(4)=fpar(1)*(nrmb+fpar(12)*nrmx)
            else
               fpar(4)=fpar(1)*nrmb
            end if
         end if
      end if
c     ||ML^{-1}r ||, actual preconditioned residual
      fpar(5) = alpha
c     leading component of the right hand side in 
c     ||ML^{-1}(b-Ax)||=||alpha e_1 - Hy||
      w(vrn+1) = alpha

c         ||ML^{-1}r||<= target     &   residual norm 
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999 
     +     .and. ipar(3).ne.3) then
c        exit
         ipar(1) = 0
c        norm ||ML^{-1}(b-Ax)|| on exit
         fpar(6) = alpha
         goto 300
      elseif (alpha0.le.fpar(4) .and. ipar(3).eq.3) then
c        exit
         ipar(1) = 0
c        norm ||(b-Ax)|| on exit
         fpar(6) = alpha0
         goto 300
      endif
c     reciprocal for V(:,1)=ML^{-1}r/alpha
      alpha = one / alpha

c     normalize w
      j=1
      beta=alpha
      call SCAL(n,beta,w(1),j)

c     flop counter
      fpar(11) = fpar(11) + n
c

c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c                         MAIN LOOP
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c     request for (1) right preconditioning
c                 (2) matrix vector multiplication
c                 (3) left preconditioning
c
c     increase number of steps until restart
 110  k = k + 1
c
c     w(:,k+1)=ML^{-1}AMR^{-1}w(:,k)
c
c     right preconditioning
      if (rp) then
c        request for right preconditioning v=MR^{-1}u
         ipar(1) = 5
c        source u=w(:,k)
         ipar(8) = k*n - n + 1
         if (lp) then
c           result v=w(:,k+1)
            ipar(9) = k*n + 1
         else
c           result z=w(:,m+2), m+1 vectors are needed for restart
            ipar(9) = idx + 1
         endif
c        return flag, indicate that we want to re-enter at the third label 30
         ipar(10) = 3
         return
      endif
c     ------------------------------------------------------------------

c     ------------------------------------------------------------------
c     request for x=Av
 30   ipar(1) = 1
c      write (6,'(A)') 'return from prec., apply mat-vec'
      if (rp) then
c        if right preconditioning has been applied, then source v is located
c        at the position of the result v from v=MR^{-1}u
         ipar(8) = ipar(9)
      else
c        source v=w(:,k)
         ipar(8) = (k-1)*n + 1
      endif
c     if left preconditioning is also desired 
      if (lp) then
c        result x=w(:,m+2), m+1 vectors are needed for restart
         ipar(9) = idx + 1
      else
c        no left prec., result x=w(:,k+1), space for the next Arnoldi vector
         ipar(9) = 1 + k*n
      endif
c     return flag, indicate that we want to re-enter at the fourth label 40
      ipar(10) = 4
      return
c     ------------------------------------------------------------------
c
c     ------------------------------------------------------------------
c     left preconditioning y=ML^{-1}x
 40   if (lp) then
c        request for left preconditioning y=ML^{-1}x
         ipar(1) = 3
c        source x is obtained from previous step x=Av
         ipar(8) = ipar(9)
c        result y is returned in y=w(:,k+1)
         ipar(9) = k*n + 1
c        return flag, indicate that we want to re-enter at the fifth label 50
         ipar(10) = 5
         return
      endif
c      write (6,'(A)') 'return from mat-vec., compute Gram-Schmidt'

c     ------------------------------------------------------------------
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Hessenberg matrix. 'p2' points to the new basis vector
c
c     increment step counter
 50   ipar(7) = ipar(7) + 1
c     space to column of H(:,k+1), the previous columns are already transformed
c     to triangular form using Givens rotations
      ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
c     reorthogonalize w(:,k+1) against its predecessors using MGS with
c     selective reorthogonalization
      call MGSRO(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),
     $     ipar(12))
c      write (6,'(A,I4)') 'Gram-Schmidt completed, ipar(12)=',ipar(12)
      if (ipar(12).lt.0) goto 200
c
c     apply previous Givens rotations and generate a new one to eliminate
c     the subdiagonal element.
c
c     position of H(:,k)
      p2 = ptr + 1
c     H(1:k,k) is transformed using of old rotations
      ii=1
      do i = 1, k-1
c        H(i:i+1,k) is rotated
         ptr = p2
         p2 = p2 + 1
         call ROT(ii,w(ptr),ii,w(p2),ii,w(vc+i),w(vs+i))
      enddo
c      write (6,'(A)') 'last column of H updated'
C     eliminate H(k:k+1,k) by Givens rotation
      call ROTG(w(p2), w(p2+1), c, s)
c     store c and s
      w(vc+k) = c
      w(vs+k) = s
c     update right hand side
c     [alpha; 0] -> [c; -s']*alpha
      p2 = vrn + k
      beta = - CONJ(s) * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = beta
c      write (6,'(A)') 'end Arnoldi'
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
c     flop counter
      fpar(11) = fpar(11) + 6*k + 2
c     new residual norm
      alpha = ABS(beta)
      fpar(5) = alpha

c     k<m:                          within a restarted sequence
c     ipar(3)>=0 & alpha<=fpar(4):  ||ML^{-1}r||<=threshold (relative and/or absolute)
c     ipar(6)<=0:                   no max. its.
c     ipar(7)<ipar(6):              number of iterations less than max. its.
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.adjust*fpar(4))
     +     .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c                         END MAIN LOOP
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Hessenberg matrix,
c     p2 points to the right-hand-side (also the solution) of the system.
c
c     ------------------------------------------------------------------
 200  ptr = hess + k * (k + 1) / 2
      p2 = vrn + k
      if (w(ptr).eq.zero) then
c
c     if the diagonal elements of the last column is zero, reduce k by 1
c     so that a smaller triangular system is solved [It should only
c     happen when the matrix is singular, and at most once!]
c
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
c     H -> Givens -> R
c     backward solve R*y=w(vrn+(1:k))
      w(p2) = w(p2) / w(ptr)
      j=1
      do i = k-1, 1, -1
         ptr = ptr - i - 1

         w(p2)=-w(p2)
         call AXPY(i,w(p2),w(ptr+1),j,w(vrn+1),j)
         w(p2)=-w(p2)

         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
c     ------------------------------------------------------------------
c
c     update approximate solution
c     
c     ------------------------------------------------------------------
c     first compute w(:,1)=Vy
c
      call SCAL(n,w(p2),w(1),j)

      do i = 1, k-1
         ptr = i*n
         p2 = p2 + 1

         call AXPY(n,w(p2),w(ptr+1),j,w(1),j)
      enddo
c     flop counter
      fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
c     ------------------------------------------------------------------
c
c     ------------------------------------------------------------------
c     replace w(:,1)=Vy by MR^{-1} w(:,1)
      if (rp) then
c        request for right preconditioning
         ipar(1) = 5
c        source w(:,1)
         ipar(8) = 1
c        result w(:,m+2)
         ipar(9) = idx + 1
c        return to sixth label 60
         ipar(10) = 6
         return
      endif
c     ------------------------------------------------------------------
c
c     ------------------------------------------------------------------
c     update approximate solution
 60   if (rp) then
         do i = 1, n
c           update solution using MR^{-1} Vy
            sol(i) = sol(i) + w(idx+i)
         end do
         if (iabs(ipar(3)).eq.3) then
            i=1
            deltax=NRM(n,w(idx),i)
         end if
      else
         do i = 1, n
c           without right preconditioning we can directly update sol using w(:,1)
            sol(i) = sol(i) + w(i)
         enddo
         if (iabs(ipar(3)).eq.3) then
            i=1
            deltax=NRM(n,w,i)
         end if
      endif
c     increase flop counter
      fpar(11) = fpar(11) + n
c     ------------------------------------------------------------------
c
c     process the complete stopping criteria
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      elseif (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = ABS(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         fpar(6) = ABS(w(vrn+k))
      elseif (ipar(3).eq.1 .or. ipar(3).eq.2) then
c        store norm of the preconditioned residual
         fpar(6) = fpar(5)
      else
         i=1
c        update right hand side for the backward error
         nrmx=NRM(n,sol,i)
         if (deltax.le.fpar(1)*nrmx) then
            fpar(4)=fpar(1)*(nrmb+fpar(12)*nrmx)
         else
            fpar(4)=fpar(1)*nrmb
         end if
c         fpar(4)=fpar(1)*(nrmb+fpar(12)*nrmx)

c        ------------------------------------------------------------------
c        compute exact residual b-Ax

c        set status flag for a request for performing y=A*x
         ipar(1) = 1
c        where does the source x start: w(:,2)
         ipar(8) = n+1
c        where do we ask for the return value y: w(:,1)
         ipar(9) = 1
c        return flag, indicate that we want to re-enter at the eighth label 80
         ipar(10) = 8
c        reset counter for a restarted GMRES sequence      
         k = 0

c        copy the initial solution to the buffer for the source x
c        w(:,2)=sol
         j=1
         call COPY(n,sol,j,w(n+1),j)
         return
c        ------------------------------------------------------------------

      endif
c
c     do we need to restart ?
c
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and.
     +     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
     +     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
c
c     termination, set error code, compute convergence rate
c
 310  if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
     +     ipar(7).gt.ipar(13)) then
         fpar(7) = LOG10(fpar(3) / fpar(6)) / IREAL(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      return
c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c     return from matrix-vector multiplication
 80   if (lp) then
c        build residual r=b-Ax, store it in w(:,2), this is the place
c        for the source for left preconditioning
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         end do
c        alpha0 = ||r||
         i=1
         alpha0 = NRM(n,w(n+1),i)

c        flop counter
         fpar(11) = fpar(11) + 3*n

c        set status flag for a request for z=ML^{-1}*r
         ipar(1) = 3
c        return flag, indicate that we want to re-enter at the second label 20
         ipar(10) = 2
      else
c        build residual r=b-Ax, rewrite it to w(:,1), this is the place where
c        left preconditioning would have returned the result
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
c        alpha0 = ||r||
         i=1
         alpha0 = NRM(n,w,i)

c        flop counter
         fpar(11) = fpar(11) + 3*n
      endif
c     ------------------------------------------------------------------

      if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
c     store norm of the exact residual
      fpar(6)=alpha0
c     adjust relation between the original and the preconditioned residual
c     only the preconditioned residual can be easliy measured during a single
c     restart sequence. For safety, add a factor 0.5
      adjust=0.5* (fpar(5)/alpha0)
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and.
     +     fpar(6).gt.fpar(4)) then
c        restart GMRES

c        increment step counter
         ipar(7) = ipar(7) + 1
         ipar(13) = ipar(13) + 1
         if (lp) then
c           call for left preconditioning to re-enter restarted GMRES
            return
         else
c           re-enter restarted GMRES without left preconditioning
            goto 20
         end if
      end if
      goto 310
      end
c-----end-of-gmres
