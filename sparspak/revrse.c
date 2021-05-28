#include <sparspak.h>

/*****************************************************************************/
/* 'revrse' flips the components of the vector 'v', the number of components */
/* of 'v' is 'nv'. This function is cited in the SPARSPAK functions 'gennd', */
/* 'gen1wd' but it is not included in the source code.                       */
/* in C, 'v' has subscripts 0,1,2,...,*nv-1.                                 */
/*****************************************************************************/
void revrse(integer *nv, integer *v)
{
   integer   dummy, 
         i=*nv/2, k;

   k=i%4;
   while (i>k)
   {     
         /* exchange entries */
         dummy=v[i-1];   v[i-1]=v[*nv-i];     v[*nv-i]  =dummy;
         dummy=v[i-2];   v[i-2]=v[*nv-i+1];   v[*nv-i+1]=dummy;
         dummy=v[i-3];   v[i-3]=v[*nv-i+2];   v[*nv-i+2]=dummy;
         dummy=v[i-4];   v[i-4]=v[*nv-i+3];   v[*nv-i+3]=dummy;
	 i-=4;
   } /* end while */
   while (i)
   {     
         /* exchange entries */
         dummy=v[i-1];	 v[i-1]=v[*nv-i];     v[*nv-i]  =dummy;
	 i--;
   } /* end while */
} /* end revrse */
