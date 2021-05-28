#include <ilupack_fortran.h>
#include <ilupackmacros.h>



      logical function STOPBIS(n,ipar,mvpi,fpar,r,delx,sx)
      implicit none
      integer n,mvpi,ipar(16)
      REALS fpar(16), sx, NRM
      FLOAT r(n), delx(n)
      external NRM
c-----------------------------------------------------------------------
c     function for determining the stopping criteria. return value of
c     true if the stopbis criteria is satisfied.
c
c     code taken from SPARSKIT of Yousef Saad.
c     adapted by Matthias Bollhoefer for the complex case
c-----------------------------------------------------------------------
      integer i

c      write (6,'(A,I3)')'stopbis, on entry ipar(1)=',ipar(1)
      i=1
      if (ipar(11) .eq. 1) then
         STOPBIS = .true.
      else
         STOPBIS = .false.
      endif
c     ipar(6): maximum number of steps
c     ipar(7): current number of steps
      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
c        maximum number of steps exceeded
         ipar(1) = -1
         STOPBIS = .true.
      endif
      if (STOPBIS) then 
c         write (6,'(A,I3)')'stopbis, on exit ipar(1)=',ipar(1)
         return
      end if
c
c     computes errors
c
c     current norm of the residual
      fpar(5) = NRM(n,r,i)
c     increment flop counter
      fpar(11) = fpar(11) + 2 * n
c     type of stopping criterion
      if (ipar(3).lt.0) then
c
c     compute the change in the solution vector
c
         fpar(6) = sx * NRM(n,delx,i)
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).lt.mvpi+mvpi+1) then
c
c     if this is the end of the first iteration, set fpar(3:4)
c
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * fpar(3) + fpar(2)
            endif
         endif
      else
c        transfer norm of the residual
         fpar(6) = fpar(5)
      endif
c
c     .. the test is struct this way so that when the value in fpar(6)
c       is not a valid number, STOPBIS is set to .true.
c
      if (fpar(6).gt.fpar(4)) then
         STOPBIS = .false.
         ipar(11) = 0
      else
         STOPBIS = .true.
         ipar(11) = 1
      endif
c
c      write (6,'(A,I3)')'stopbis, on exit ipar(1)=',ipar(1)
      return
      end
c-----end-of-stopbis
