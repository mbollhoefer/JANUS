#include <ilupack_fortran.h>
      double complex function zdotc2(n,zx,incx)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
c     changes done by Matthias Bollhoefer, October 2003 
c     copy of zdotc when x=y. This version avoids conflicts
c     with the optimization when alias of x and y are used
c     
      implicit none
      double complex zx(1),ztemp
      integer i,incx,ix,n
      ztemp = (0.0d0,0.0d0)
      zdotc2 = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zx(ix)
        ix = ix + incx
   10 continue
      zdotc2 = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zx(i)
   30 continue
      zdotc2 = ztemp
      return
      end
