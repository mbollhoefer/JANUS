#include <ilupack_fortran.h>
#include <ilupackmacros.h>



      subroutine MGSRO(full,lda,n,m,ind,ops,vec,hh,ierr)
      implicit none
      logical full
      integer lda,m,n,ind,ierr
      REALS ops
      FLOAT  hh(m),vec(lda,m)
c-----------------------------------------------------------------------
c     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
c               Orthogonalization
c     The ind'th vector of VEC is orthogonalized against the rest of
c     the vectors.
c
c     The test for performing re-orthogonalization is performed for
c     each indivadual vectors. If the cosine between the two vectors
c     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
c     performed. The norm of the 'new' vector is kept in variable NRM0,
c     and updated after operating with each vector.
c
c     full   -- .true. if it is necessary to orthogonalize the ind'th
c               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
c               .false. only orthogonalize againt vec(:,1:ind-1)
c     lda    -- the leading dimension of VEC
c     n      -- length of the vector in VEC
c     m      -- number of vectors can be stored in VEC
c     ind    -- index to the vector to be changed
c     ops    -- operation counts
c     vec    -- vector of LDA X M storing the vectors
c     hh     -- coefficient of the orthogonalization
c     ierr   -- error code
c               0 : successful return
c               -1: zero input vector
c               -2: input vector contains abnormal numbers
c               -3: input vector is a linear combination of others
c
c     External routines used: FLOAT DISTDOT
c
c     code taken from SPARSKIT of Yousef Saad.
c     adapted by Matthias Bollhoefer for the complex case
c-----------------------------------------------------------------------
      integer i,k,j,jj
      REALS  nrm0, nrm1, thr, zero, one, reorth, SQRTFNAME, NRM
      FLOAT fct
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_
      parameter (zero=0.0e0, one=1.0e0, reorth=0.98e0)
#else
      parameter (zero=0.0D0, one=1.0D0, reorth=0.98D0)
#endif
      external DISTDOT, DISTDOT2,AXPY,NRM
      FLOAT DISTDOT, DISTDOT2
      intrinsic SQRTFNAME
c
c     compute the norm of the input vector
c
      j=1
      jj=1
      nrm0 = CTOD(DISTDOT2(n,vec(1,ind),j))
c      write (6,'(A,1P,E12.4)') 'start Gram-Schmidt nrm0=',nrm0
      ops = ops + n + n
      thr = nrm0 * reorth
      if (nrm0.le.zero) then
         ierr = - 1
         return
      else if (nrm0.gt.zero .and. one/nrm0.gt.zero) then
         ierr = 0
      else
         ierr = -2
         return
      endif
c
c     Modified Gram-Schmidt loop
c
      if (full) then
c         write (6,'(A,I4,A,I4)')
c     +        '1st Gram-Schmidt loop',ind+1,',...,',m
         do 40 i = ind+1, m
c            write (6,'(A,I4)') ' Gram-Schmidt step',i
            fct = DISTDOT(n,vec(1,i),j,vec(1,ind),jj)
            hh(i) = fct
            fct=-fct
            call AXPY(n,fct,vec(1,i),j,vec(1,ind),jj)
            fct=-fct
            ops = ops + 4 * n + 2
            if (CTOD(fct*CONJ(fct)).gt.thr) then
               fct = DISTDOT(n,vec(1,i),j,vec(1,ind),jj)
               hh(i) = hh(i) + fct
               fct=-fct
               call AXPY(n,fct,vec(1,i),j,vec(1,ind),jj)
               fct=-fct
               ops = ops + 4*n + 1
            endif
            nrm0 = nrm0 - CTOD(hh(i) * CONJ(hh(i)))
            if (nrm0.lt.zero) nrm0 = zero
            thr = nrm0 * reorth
 40      continue
      endif
c
c      write (6,'(A,I4,A,I4)') 
c     +     '2nd Gram-Schmidt loop',1,',...,',ind-1
      do 70 i = 1, ind-1
c         write (6,'(A,I4)') '2nd Gram-Schmidt step',i
c         write (6,'(1P,E12.4)') dznrm2(n,vec(1,ind),j)
c         write (6,'(1P,E12.4)') dznrm2(n,vec(1,i),j)
         fct = DISTDOT(n,vec(1,i),j,vec(1,ind),jj)
c         write (6,'(A,1P,E12.4)') 'fct=',fct
         hh(i) = fct
         fct=-fct
         call AXPY(n,fct,vec(1,i),j,vec(1,ind),jj)
c         write (6,'(A)') 'axpy completed'
         fct=-fct
         ops = ops + 4 * n + 2
         if (CTOD(fct*CONJ(fct)).gt.thr) then
            fct = DISTDOT(n,vec(1,i),j,vec(1,ind),jj)
c            write (6,'(A,1P,E12.4)') 'fct=',fct
            hh(i) = hh(i) + fct
            fct=-fct
            call AXPY(n,fct,vec(1,i),j,vec(1,ind),jj)
c            write (6,'(A)') 'axpy completed'
            fct=-fct
            ops = ops + 4*n + 1
         endif
         nrm0 = nrm0 - CTOD(hh(i) * CONJ(hh(i)))
         if (nrm0.lt.zero) nrm0 = zero
         thr = nrm0 * reorth
 70   continue
c
c     test the resulting vector
c
      nrm1 = SQRT(CTOD(DISTDOT2(n,vec(1,ind),j)))
c      write (6,'(A,1P,E12.4)') ' norm resulting vector',nrm1
      ops = ops + n + n
 75   hh(ind) = nrm1
      if (nrm1.le.zero) then
         ierr = -3
         return
      endif
c
c     scale the resulting vector
c
      fct = one / nrm1
      call SCAL(n,fct,vec(1,ind),j)
      ops = ops + n + 1
c
c     normal return
c
      ierr = 0
      return
c     end subroutine mgsro
      end
