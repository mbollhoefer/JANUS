#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <janus.h>

#include <ilupackmacros.h>

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout


//#define PRINT_INFO




/* compute symmetric weighted matching, i.e. break cycles into
   1x1 and 2x2 cycles

   A     is a GENERAL but (almost) symmetrically structured matrix
   p     is a vector of length A.nc carrying the old and the new
         permutation
   ibuff is a buffer of length A.nc
   dbuff is a buffer of length 3*A.nc

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        July/November 2003, 2005. ILUPACK V2.1

    Notice:

	Copyright (c) 2005 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
*/
integer GNLSYMMWM(CSRMAT A, integer *p, integer *ibuff, FLOAT *dbuff)
{
  integer i,j,k,l=1,n=A.nc, m, flag, next;
  REALS weight, weight2, weighte, weighto, val, val2, 
        *rbuff=(REALS *)(dbuff+A.nc),*cbuff=(REALS *)(dbuff+2*A.nc);
  FLOAT aij, aji;
  
  for (i=0; i<n; i++) {
      cbuff[i]=0.0;
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
      dbuff[i]=0.0;
#else
      dbuff[i].r=0.0;
      dbuff[i].i=0.0;
#endif
      ibuff[i]=0;
  }

  // compute 1-norm for the off-diagonal entries by columns and by rows and 
  // extract diagonal entries
  for (i=0; i<n; i++) {
      j=A.ia[i]-1;
      k=A.ia[i+1]-1-j;
      l=1;
      rbuff[i]=ASUM(&k,A.a+j,&l);
      k+=j;
      for (; j<k; j++) {
	  l=A.ja[j]-1;
	  if (l!=i)
	     cbuff[l]+=FABS(A.a[j]);
	  else {
	     rbuff[i]-=FABS(A.a[j]);
	     dbuff[i]=A.a[j];
	  }
      } // end for j
  } // end for i

  
  // break cycles into 1x1 and 2x2 cycles
  // first element of a cycle
  l=0;
  while (l<n) {

        flag=-1;
	// move to the next undiscovered cycle
	while (flag) 
	      if (l>=n || ibuff[l]==0)
		 flag=0;
	      else
		 l++;

	// mark entries of cycle k
	if (l<n) {


#ifdef PRINT_INFO
	   printf("%8d",l+1);
#endif

	   // mark element l as member of this cycle
	   ibuff[l]=-1;
	   // successor in this cycle
	   i=p[l]-1;
	   // length of the cycle
	   j=1;
	   while (i!=l) {


#ifdef PRINT_INFO
	         printf("%8d",i+1);
#endif


		 // mark element i as member of this cycle
	         ibuff[i]=-1;
		 i=p[i]-1;
		 j++;
	   } // end while


#ifdef PRINT_INFO
	   printf("\n");
	   fflush(stdout);
#endif

	   // even cycle of length greater than 2
	   if (j>2 && j%2==0) {
	      /* break even cycle into a product of subsequent 2-cycles.
		 There are two isomorphic possibilities to break up the cycle.
		 The first sequence starts with a p[l], the second one with p[p[l]-1].
		 We break the cycle such that the more diagonal dominant sequence of
		 2-cycles is taken
	      */
	      weighte=0.0;
	      weighto=0.0;
	      flag=0;
	      i=l;
	      next=p[i]-1;
	      k=0;
	      while (k<j) {
		    // find off-diagonal entries A(i,next), A(next,i)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    aij=0.0;
		    aji=0.0;
#else
		    aij.r=aij.i=0.0;
		    aji.r=aji.i=0.0;
#endif
		    // scan row i for A(i,next)
		    for (m=A.ia[i]-1; m<A.ia[i+1]-1; m++)
		        if (A.ja[m]-1==next)
			   aij=A.a[m];
		    // scan row 'next' for A(next,i)
		    for (m=A.ia[next]-1; m<A.ia[next+1]-1; m++)
		        if (A.ja[m]-1==i)
			   aji=A.a[m];

		    // compute maximum off-diagonal 1-norm of the associated rows
		    // remove |aij| and |aji|
		    val =FABS(aij);
		    val2=FABS(aji);
		    weight =MAX(rbuff[i]-val, rbuff[next]-val2);
		    weight2=MAX(cbuff[i]-val2,cbuff[next]-val);
		    weight=MAX(weight,weight2);
		    if (weight<0.0)
		       weight=0.0;
		    /* To measure the block diagonal dominance we use
		       ||/aii aij\^{-1}||          ||/ ajj -aij\||      weight
		       |||       |     ||*weight = |||         |||* ---------------
 		       ||\aji ajj/     ||          ||\-aji  aii/|| |aii*ajj-aij*aji|  

                         MAX(|ajj|+|-aij|,|-aji|+|aii|)*weight
                       = -------------------------------------
                                   |aii*ajj-aij*aji|  

                         (MAX(|ajj|,|aii|)+|aij|)*weight
                       = -------------------------------
                                |aii*ajj-aij*aji|
		    */
		    weight*=sqrt(MAX(FABS(dbuff[i])+val2,FABS(dbuff[next])+val)
                                *MAX(FABS(dbuff[i])+val, FABS(dbuff[next])+val2));


		    // build determinant aii*ajj-aij*aji 
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		    // compute determinant
		    aij=dbuff[i]*dbuff[next]-aij*aji;
#else
		    // aij*aji
		    val=aij.r;
		    aij.r=aij.r*aji.r-aij.i*aji.i;
		    aij.i=val  *aji.i+aij.i*aji.r;

		    // compute determinant
		    aij.r=dbuff[i].r*dbuff[next].r-dbuff[i].i*dbuff[next].i-aij.r;
		    aij.i=dbuff[i].r*dbuff[next].i+dbuff[i].i*dbuff[next].r-aij.i;
#endif

		    // catch the case that the 2x2 block is singular
		    weight/=FABS(aij)+1.0e-30;

		    if (flag) {
		       weighto+=weight;


#ifdef PRINT_INFO
		       printf("o%8.1le,%8.1le\n",weight,weighto);
#endif


		       flag=0;
		    }
		    else {
		       weighte+=weight;


#ifdef PRINT_INFO
		       printf("e%8.1le,%8.1le\n",weight,weighte);
#endif


		       flag=-1;
		    }
		    i=next;
		    next=p[i]-1;
		    k++;
	      } // end while (k<j)

	      // we prefer the sequence of 2-cycles that is more block diagonal dominant
	      if (weighte<=weighto)
		 i=l;
	      else
		 i=p[l]-1;
	      next=p[i]-1;
	      k=0;
	      while (k<j) {


#ifdef PRINT_INFO
	            printf("%8d%7d,",i+1,next+1);
#endif


		    // store i
		    m=i;
		    // advance i by two steps
		    i=p[next]-1;
		    // 2-cycle (i,next)
		    p[next]=m+1;
		    next=p[i]-1;
		    k+=2;
	      } // end while (k<j)


#ifdef PRINT_INFO
	      printf("\n\n");
	      fflush(stdout);
#endif


	   } // end if (j>2 && j%2==0)


	   // odd cycle of length greater than 2
	   else if (j>2) {
	      // find 1x1 cycle that is most diagonal dominant
	      i=l;
	      next=p[i]-1;
	      weight=1e30;
	      m=-1;
	      k=0;
	      while (k<j) {
	            val=FABS(dbuff[i]);
		    if (val!=0.0) {
		       val=MAX(rbuff[i],cbuff[i])/val;


#ifdef PRINT_INFO
		       printf("%8.1le\n",val);
#endif


		       if (val<weight) {
			  weight=val;
			  m=i;
		       }
		    }
		    i=next;
		    next=p[i]-1;
		    k++;
	      } // end while

	      // did we find at least one nonzero diagonal entry?
	      if (m>=0) {
		 // we take the most diagonal dominant diagonal entry as 1x1 cycle
		 // break the cycle into a sequence of 2-cycles behind the singleton
		 i=p[m]-1;
		 // singleton
		 p[m]=m+1;


#ifdef PRINT_INFO
		 printf("%7d,",m+1);
#endif


		 next=p[i]-1;
		 k=1;
		 while (k<j) {


#ifdef PRINT_INFO
		       printf("%8d%7d,",i+1,next+1);
#endif


		       // store i
		       m=i;
		       // advance i by two steps
		       i=p[next]-1;
		       // 2-cycle p[i] and p[next]
		       p[next]=m+1;
		       next=p[i]-1;
		       k+=2;
		 } // end while


#ifdef PRINT_INFO
		 printf("\n\n");
		 fflush(stdout);
#endif


		 
	      } // end if (m>=0)
	      else { // all diagonal entries are zero 
		 /* only check two subsequent cycle sequences similar to the even case
		    break cycle such that the more diagonal dominant sequence is taken.
		    We ignore the final diagonal entry (which is zero in any case)
		 */
		 weighte=0.0;
		 weighto=0.0;
		 flag=0;
		 i=l;
		 next=p[i]-1;
		 k=1;
		 while (k<j) {
		       // find off-diagonal entries A(i,next), A(next,i)
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		       aij=0.0;
		       aji=0.0;
#else
		       aij.r=aij.i=0.0;
		       aji.r=aji.i=0.0;
#endif
		       // scan row i for A(i,next)
		       for (m=A.ia[i]-1; m<A.ia[i+1]-1; m++)
			   if (A.ja[m]-1==next)
			      aij=A.a[m];
		       // scan row 'next' for A(next,i)
		       for (m=A.ia[next]-1; m<A.ia[next+1]-1; m++)
		           if (A.ja[m]-1==i)
			      aij=A.a[m];

		       // compute maximum off-diagonal 1-norm of the associated rows
		       // remove |aij| and |aji| which have the same absolute value
		       val =FABS(aij);
		       val2=FABS(aji);
		       weight =MAX(rbuff[i]-val, rbuff[next]-val2);
		       weight2=MAX(cbuff[i]-val2,cbuff[next]-val);
		       weight=MAX(weight,weight2);
		       if (weight<0.0)
			  weight=0.0;
		       /* To measure the block diagonal dominance we use
			  ||/aii aij\^{-1}||          ||/ ajj -aij\||      weight
			  |||       |     ||*weight = |||         |||* ---------------
			  ||\aji ajj/     ||          ||\-aji  aii/|| |aii*ajj-aij*aji|  

                             MAX(|ajj|+|-aij|,|-aji|+|aii|)*weight
			  =  -------------------------------------
                                       |aii*ajj-aij*aji|  

                            (MAX(|ajj|,|aii|)+|aij|)*weight
			  = -------------------------------
                                   |aii*ajj-aij*aji|
		       */
		       weight*=sqrt(MAX(FABS(dbuff[i])+val2,FABS(dbuff[next])+val)
				   *MAX(FABS(dbuff[i])+val, FABS(dbuff[next])+val2));

		       
		       // build determinant aii*ajj-aij*aji 
		       // To do this form product A(p[i],p[next])*A(p[next],p[i])
#if defined _SINGLE_REAL_ || defined _DOUBLE_REAL_
		       // compute determinant
		       aij=dbuff[i]*dbuff[next]-aij*aji;
#else
		       // aij*aji
		       val=aij.r;
		       aij.r=aij.r*aji.r-aij.i*aji.i;
		       aij.i=val  *aji.i+aij.i*aji.r;

		       // compute determinant
		       aij.r=dbuff[i].r*dbuff[next].r-dbuff[i].i*dbuff[next].i-aij.r;
		       aij.i=dbuff[i].r*dbuff[next].i+dbuff[i].i*dbuff[next].r-aij.i;
#endif

		       // catch the case that the 2x2 block is singular
		       weight/=FABS(aij)+1.0e-30;

		       if (flag) {
			  weighto+=weight;


#ifdef PRINT_INFO
		       printf("o%8.1le,%8.1le\n",weight,weighto);
#endif


			  flag=0;
		       }
		       else {
			  weighte+=weight;


#ifdef PRINT_INFO
		       printf("e%8.1le,%8.1le\n",weight,weighte);
#endif


			  flag=-1;
		       }
		       i=next;
		       next=p[i]-1;
		       k++;
		 } // end while (k<j)

		 // we choose the cycle that is more block diagonal dominant
		 if (weighte<=weighto)
		    i=l;
		 else
		    i=p[l]-1;
		 next=p[i]-1;
		 k=1;
		 while (k<j) {


#ifdef PRINT_INFO
		       printf("%8d%7d,",i+1,next+1);
#endif


		       // store i
		       m=i;
		       // advance i by two steps
		       i=p[next]-1;
		       // 2-cycle (i,next)
		       p[next]=m+1;
		       next=p[i]-1;
		       k+=2;
		 } // end while
		 // finally set singleton
		 p[i]=i+1;


#ifdef PRINT_INFO
		 printf("%7d,\n\n",i+1);
		 fflush(stdout);
#endif


	      } // end if-else (m>=0)
	   } // end if-else-if (j>2)
	} // end if  (l<n)
  } // end while (l<n)
  
  // finally rearrange the permutation such that the main weight is moved
  // to the tridiagonal part of A
  for (i=0; i<n; i++)
      ibuff[i]=p[i];

  j=0;
  for (i=0; i<n; i++) {
      k=ibuff[i];
      if (k-1==i)
	 p[j++]=k;
  } // end for i
  l=j;
  for (i=0; i<n; i++) {
      k=ibuff[i];
      if (k-1!=i) {
	 p[j++]=k;
	 p[j++]=i+1;
	 ibuff[k-1]=k;
      }
  } // end for i
  return (l);
} // end swm

