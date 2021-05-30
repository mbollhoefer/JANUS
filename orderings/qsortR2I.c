#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <janus.h>

#include <ilupackmacros.h>


void QSORTR2I(REALS *wa, integer *cor1, integer *cor2, integer left, integer right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into decreasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   integer i, last;

   if (left >= right)  return;

   SWAPM(wa, left, (left+right)/2);
   swapj(cor1, left, (left+right)/2);
   swapj(cor2, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] > wa[left]) {
	 SWAPM(wa, ++last, i);
	 swapj(cor1, last, i);
	 swapj(cor2, last, i);
      }
   }
   SWAPM(wa, left, last);
   swapj(cor1, left, last);
   swapj(cor2, left, last);
   QSORTR2I(wa, cor1, cor2, left, last-1);
   QSORTR2I(wa, cor1, cor2, last+1, right);
}
