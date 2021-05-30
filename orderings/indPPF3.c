#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <janus.h>


#include <ilupackmacros.h>

#define  ALPHA  0.00001





void SPDWDIAG(CSRMAT, integer *, REALS, integer *, REALS *, REALS *, REALS *, integer *);


integer PPPERM(CSRMAT mat, integer bsize, integer *Qord, integer *nnod, 
	   REALS tol, FLOAT *dbuff, integer *ibuff) {
/*--------------------------------------------------------------------- 
| algorithm for symmetric  block selection - 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     mat  =  matrix in compressed sparse row format
|     
|     tol    =  a tolerance for excluding a row from B block 
|
|     bsize not used here - it is used in arms2.. 
|
|     Output parameters:
|     ------------------ 
|     Qord   = column/row permutation array.  Column number j will become 
|	       row number Qord[j] in permuted matrix.
|             [destination lists] - i.e., old to new arrays= 
|     
|     nnod   = number of elements in the B-block 
|     
|----------------------------------------------------------------------
| original unsymmetric version written by Yousef Saad, 2003
| adaption to the symmetric case by Matthias Bollhoefer, November 2003
|---------------------------------------------------------------------*/ 
/*--------------------   local variables   */
  integer *jcor;
  integer i, j,  n=mat.nr, nnzrow, count, numnode;
  integer ii, jj, k, exclu, *row; 
  REALS wmax,wsum,s,t, mx;
  FLOAT *mrow;

/*-------------------- prototypes */
  integer check_perm(integer, integer *) ;   /* debug only */
/*-----------------------------------------------------------------------*/  
  
   for (j=0; j<n; j++) {
     Qord[j] = -1; 
   } // end for j

   jcor = ibuff;     // n   spaces required 

/*--------------------------------------------------
|  get coordinate format
|-------------------------------------------------*/
  numnode = 0;
  count = 0;
/*-------------------- wDiag selects candidate entries in sorted oder */
  /* use a weaker drop tolerance when preselecting rows/cols with a certain
     'diagonal dominance' weight */
  SPDWDIAG(mat, jcor, tol, &count, &wmax,&wsum, (REALS *)dbuff, ibuff+n); 
/*-------------------- add entries one by one to diagnl */
  for (i = 0; i<count; i++){ 
    jj=jcor[i]; 
    Qord[jj] = numnode;
    numnode++; 
  }
  /*-------------------- number of B nodes */
  *nnod = numnode; 
/*--------------------------------------------------
|    end-main-loop - complete permutation arrays
|-------------------------------------------------*/
  for (j=0; j<n; j++)
    if (Qord[j] < 0) 
      Qord[j] = numnode++;

  if (numnode != n) {
    printf("  ** counting error - type 1 \n"); return 1; }   

/*--------------------------------------------------
|  clean up before returning
|-------------------------------------------------*/
  return 0;
}
/*---------------------------------------------------------------------
|-----end-of-indsetPP--------------------------------------------------
|--------------------------------------------------------------------*/

void SPDWDIAG(CSRMAT mat, integer *jcor, REALS tol, integer *count, 
	   REALS *wmax, REALS *wsum, REALS *dbuff, integer *ibuff) 
{
/*---------------------------------------------------------------------
|     defines weights based on diagonal dominance relative to 
|     excluded rows and columns
|    
|     fir indPO2.. [upper triangular reduction] 
|--------------------------------------------------------------------*/
  integer i, k, n=mat.nr, len, col, jmax, countL;
  integer *nz, *nzT, *jcol, *imax; 
  REALS *rownorm, *colnorm, *weight, s,t, tmax, *smax;
  FLOAT *mrow;

/*--------------------begin */

  rownorm=dbuff;     // n   spaces required, only temporarily needed
  colnorm=dbuff+n;   // n   spaces required   totally 3*n REALS spaces
  smax   =dbuff+2*n; // n   spaces required   
  weight =rownorm;   // alias on rownorm, rownorm is only needed temporarily
 
  nz     =jcor;      // alias on jcor, only temporarily needed
  nzT    =ibuff;     // n   spaces required    totally n integer    spaces

  len = 0; 
  for (i=0; i<n; i++) {
    rownorm[i]=0.0;
    colnorm[i]=0.0;
    smax[i]   =0.0;

    nz[i]  =0; 
    nzT[i] =0; 
  }
/*-------------------- compute row + column norms and nzT */
  for (i=0; i<n; i++) {
      jcol=mat.ja+mat.ia[i]-1;
      mrow=mat.a +mat.ia[i]-1;
      // absolute value of the diagonal entry
      tmax=0.0;
      for (k =0; k<mat.ia[i+1]-mat.ia[i]; k++) {
	  col = jcol[k]-1; 
	  t = FABS(mrow[k]); 
	  if (t!=0.0) {
	     rownorm[i]  +=t; 
	     nz[i]++;
	     if (col==i) {
	        tmax=t;
		smax[col]=t;
	     }
	     else { // do not count the diagonal entry twice
	        colnorm[col]+=t;
	        nzT[col]++; 
	     } // end if-else
	  } // end if t!=0
      } // end for k
  } // end for i

  // Since only half of the matrix ist stored, we have to complete the data
  *wmax = 0.0; 
  *wsum = 0.0; 
  for (col=0; col<n; col++) {
      nzT[col]+=nz[col];
      colnorm[col]+=rownorm[col];

      // ratio of diagonal dominance 
      t = smax[col]/colnorm[col]; 
      // rownorm[col] no longer needed, reuse it as weight[col]
      weight[col]=t;
      // bookmark the largest diagonal dominant value
      if (*wmax<t)  *wmax=t; 
      // average of diagonal dominance
      *wsum += t; 
      
      // nz[col] no longer needed, reinitialize jcor[col]=nz[col]=0
      jcor[col]=0;
  } // end for col
  *wsum/=(REALS)n;

  // remove the rows/columns with poorest diagonal dominance
  for (col=0; col<n; col++) {
      t=weight[col];
      if (2.0*t<tol**wmax && t<*wsum && t<=0.51) {
	 jcor[col]=-1; 
      }
  } // end for col

  // downdate diagonal dominance information
  for (i=0; i<n; i++) {
      jcol=mat.ja+mat.ia[i]-1;
      mrow=mat.a +mat.ia[i]-1;
      for (k =0; k<mat.ia[i+1]-mat.ia[i]; k++) {
	  col = jcol[k]-1; 
	  t = FABS(mrow[k]); 
	  // remove t=|a_{i,col}|
	  if (jcor[i]<0) 
	     colnorm[col]-=t;
	  if (jcor[col]<0) 
	     colnorm[i]-=t;
      } // end for k
  } // end for i

  // recompute updated weights
  for (col=0; col<n; col++) {
      // ratio of diagonal dominance 
      if (jcor[col]==0)
	 t=smax[col]/colnorm[col]; 
      else
	 t=-1.0;
      weight[col]=t;
      jcor[col]=0;
  } // end for col


/*-------------------- now select according to tol and tighter bound*/
  countL=0; 
  for (i=0; i<n; i++) {
      t =weight[i];
      // tighter bound 't<tol**wmax' instead of '2.0t<tol**wmax'
      if (t<tol**wmax && t<*wsum && t<=0.51) continue;
      weight[countL]=t*t/((REALS)(nzT[i]-1)*(nzT[i]-1)+1);
      jcor[countL]=i; 
      // now use nzT as duplicate copy of jcor
      nzT[i]=i;
      countL++;
  } // end for i
  /*-------------------- sort them  */
  QSORTR2I(weight, nzT, jcor, 0, countL-1);
  *count = countL;
}
/*---------------------------------------------------------------------
|---- end of wDiag ----------------------------------------------------
|--------------------------------------------------------------------*/
