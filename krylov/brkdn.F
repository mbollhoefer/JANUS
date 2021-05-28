#include <ilupack_fortran.h>
#include <ilupackmacros.h>



      logical function BRKDN(alpha, ipar)
      implicit none
      integer ipar(16)
      REALS alpha, beta, zero, one
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_
      parameter(one=1.0e0, zero=0.0e0)
#else
      parameter(one=1.0D0, zero=0.0D0)
#endif
c-----------------------------------------------------------------------
c     test whether alpha is zero or an abnormal number, if yes,
c     this routine will return .true.
c
c     If alpha == 0, ipar(1) = -3,
c     if alpha is an abnormal number, ipar(1) = -9.
c-----------------------------------------------------------------------
      BRKDN = .false.
      if (alpha.gt.zero) then
         beta = one / alpha
         if (.not. beta.gt.zero) then
            BRKDN = .true.
            ipar(1) = -9
         endif
      else if (alpha.lt.zero) then
         beta = one / alpha
         if (.not. beta.lt.zero) then
            BRKDN = .true.
            ipar(1) = -9
         endif
      else if (alpha.eq.zero) then
         BRKDN = .true.
         ipar(1) = -3
      else
         BRKDN = .true.
         ipar(1) = -9
      endif
      return
      end
c-----end-of-brkdn
c-----------------------------------------------------------------------
