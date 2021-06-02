#include <ilupack_fortran.h>
      subroutine  crot (n,zx,incx,zy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
c     changes done by Matthias Bollhoefer, October 2003 
c     this routine is a complex version of srot. It closes
c     a gap in the blas library
c
      real c
      complex   zx(*),zy(*),ztemp,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp  = c*zx(ix) + s       *zy(iy)
        zy(iy) = c*zy(iy) - conjg(s)*zx(ix)
        zx(ix) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = c*zx(i) + s       *zy(i)
        zy(i) = c*zy(i) - conjg(s)*zx(i)
        zx(i) = ztemp
   30 continue
      return
      end
