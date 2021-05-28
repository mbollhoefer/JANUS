#include <ilupack_fortran.h>
#include <ilupackmacros.h>



      subroutine TIDYCG(n,ipar,fpar,sol,delx)
      implicit none
      integer i,n,ipar(16)
      REALS fpar(16)
      FLOAT sol(n),delx(n)
c-----------------------------------------------------------------------
c     Some common operations required before terminating the CG routines
c
c     code taken from SPARSKIT of Yousef Saad.
c     adapted by Matthias Bollhoefer for the complex case
c-----------------------------------------------------------------------
      REALS LOGTENFNAME, zero
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_
      parameter(zero=0.0e0)
#else
      parameter(zero=0.0D0)
#endif
c
c      write (6,'(A,I3)')'tidycg, on entry ipar(1)=',ipar(1)

      if (ipar(12).ne.0) then
         ipar(1) = ipar(12) 
      else if (ipar(1).gt.0) then
         if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or.
     +            fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(3).ne.3) then
            ipar(1) = -10
         else 
            ipar(1) = 0
         endif
      endif
      if (fpar(3).gt.zero .and. fpar(6).gt.zero .and.
     +     ipar(7).gt.ipar(13)) then
         fpar(7) = LOG10(fpar(3) / fpar(6)) / IREAL(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      do i = 1, n
         sol(i) = sol(i) + delx(i)
      enddo

c      write (6,'(A,I3)')'tidycg, on exit ipar(1)=',ipar(1)

      return
      end
c-----end-of-tidycg
