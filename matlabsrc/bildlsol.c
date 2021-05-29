/* $Id: bildlsol.c 4169 2018-04-26 16:34:10Z bolle $ */
/* ========================================================================== */
/* === bildlsol mexFunction ============================================== */
/* ========================================================================== */

/*
    Usage:

    Given a Hermitian block factorization BL BiD^{-1} BL^H, solve a linear system

    Example:

    z=bildlsol(BL,BiD,pivots,y)


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	September 20, 2016. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2016 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/
*/

/* ========================================================================== */
/* === Include files and prototypes ========================================= */
/* ========================================================================== */

#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <janus.h>
#include <ilupackmacros.h>
#include <lapack.h>

#define MAX(A,B) (((A)>=(B))?(A):(B))
#define MIN(A,B) (((A)>=(B))?(B):(A))
/* #define PRINT_CHECK */
/* #define PRINT_INFO  */

/* ========================================================================== */
/* === mexFunction ========================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
    mxArray            *Block, *BlockJ, *BlockI, *BlockE, *BlockD, 
                       *BL_input, *BiD_input,  *y_input, *x_output;

    mwIndex            *ja,  /* row indices of input matrix y (if sparse) */
                       *ia;  /* column pointers of input matrix y (if sparse) */
    
    integer            i,j,k,l,m,n,*p,isreal,nblocks,*pivots,npivots,invert_blocks, spd=0;

    double             *rhssolr, *pr,*pi, rval;
    doublecomplex      *rhssolc, *pz, *pD, *pE, val;
    DSparseBlockMatrix BLR,BiDR;
    ZSparseBlockMatrix BLC,BiDC;
    static int cnt_bildl=0;
    
    
    if (nrhs!=4)
       mexErrMsgTxt("Four input arguments required.");
    else if (nlhs!=1)
       mexErrMsgTxt("wrong number of output arguments.");


    /* The first input must be a cell array */
    BL_input=(mxArray *)prhs[0];
    if (!mxIsCell(BL_input)) {
       mexErrMsgTxt("bildlsol: first input matrix must be in cell format.") ;
    }
    /* get size of input cell array BL */
    nblocks=mxGetM(BL_input)*mxGetN(BL_input);

    /* The second input must be a cell array */
    BiD_input=(mxArray *)prhs[1];
    if (!mxIsCell(BiD_input)) {
       mexErrMsgTxt("bildlsol: second input matrix must be in cell format.") ;
    }

    /* The fourth input must be the r.h.s. */
    y_input=(mxArray *)prhs[3];
    if (!mxIsNumeric(y_input))
       mexErrMsgTxt("bildlsol: fourth input argument must be numeric.") ;
    n=mxGetM(y_input);
    m=mxGetN(y_input);    
    
    /* The third input must be invert_blocks */
    y_input=(mxArray *)prhs[2];
    if (!mxIsNumeric(y_input))
       mexErrMsgTxt("bildlsol: third input argument must be numeric.");
    npivots=mxGetM(y_input)*mxGetN(y_input);
    if (npivots==0) {
       pivots=NULL;
       invert_blocks=1;
    }
    else {
       pivots=(integer *)malloc((size_t)npivots*sizeof(integer));
       pr=mxGetPr(y_input);
       for (i=0; i<npivots; i++)
	   pivots[i]=pr[i];
       invert_blocks=0;
       /* LL^T=A, Cholesky case */
       if (pivots[0]==-(n+1))
	  spd=1;
    }
    /* mexPrintf("npivots=%ld,spd=%ld,invert_blocks=%ld\n",npivots,spd,invert_blocks); */
    
    /* The fourth input must be the r.h.s. */
    y_input=(mxArray *)prhs[3];

    
#ifdef PRINT_INFO
    mexPrintf("bildlsol: npivots=%ld,invert_blocks=%ld",npivots,invert_blocks);fflush(stdout);
    mexPrintf("bildlsol: check real/non-real\n");fflush(stdout);
#endif
    /* check whether the real or complex case has to be applied */
    isreal=-1;
    if (mxIsComplex(y_input))
       isreal=0;
    else {
       /* check cell arrays for non-real blocks */
       for (i=0; i<nblocks; i++) {
	   Block=mxGetCell(BiD_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BiD{i} must be a structure.");
	   BlockD=mxGetField(Block,0,"D");
	   if (BlockD==NULL)
	      mexErrMsgTxt("Field BiD{i}.D does not exist.");
	   if (!mxIsNumeric(BlockD))
	      mexErrMsgTxt("Field BiD{i}.D must be numerical.");
	   if (mxIsComplex(BlockD)) {
	      isreal=0;
	      break;
	   }

	   Block=mxGetCell(BL_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BL{i} must be a structure.");
	   BlockE=mxGetField(Block,0,"L");
	   if (BlockE==NULL)
	      mexErrMsgTxt("Field BL{i}.L does not exist.");
	   if (!mxIsNumeric(BlockE))
	      mexErrMsgTxt("Field BL{i}.L must be numerical.");
	   if (mxIsComplex(BlockE)) {
	      isreal=0;
	      break;
	   }
	   
       } /* end for i */
    } /* end else */
#ifdef PRINT_INFO
    if (isreal)
       mexPrintf("bildlsol: real case\n");
    else
       mexPrintf("bildlsol: non-real case\n");
    fflush(stdout);
#endif
    
    if (isreal) {
       /* transform data structures and right hand sides */
       BiDR.nblocks=nblocks;
       BiDR.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
       BiDR.nblockrow=NULL;
       BiDR.colind   =NULL;
       BiDR.rowind   =NULL;
       BiDR.valD     =(double **) malloc((size_t)nblocks*sizeof(double *));
       BiDR.valE=NULL;

       BLR.nblocks=nblocks;
       BLR.nblockcol=NULL;
       BLR.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
       BLR.colind   =NULL;
       BLR.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
       BLR.valE     =(double **) malloc((size_t)nblocks*sizeof(double *));
       BLR.valD=NULL;


       
       /* map input cell array to struct */
       for (i=0; i<nblocks; i++) {
	   Block=mxGetCell(BiD_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BiD{i} must be a structure.");
	   
	   BlockJ=mxGetField(Block,0,"J");
	   if (BlockJ==NULL)
	      mexErrMsgTxt("Field BiD{i}.J does not exist.");
	   if (!mxIsNumeric(BlockJ))
	      mexErrMsgTxt("Field BiD{i}.J must be numerical.");
	   k=mxGetN(BlockJ)*mxGetM(BlockJ);
	   BiDR.nblockcol[i]=k;

	   BlockD=mxGetField(Block,0,"D");
	   if (BlockD==NULL)
	      mexErrMsgTxt("Field BiD{i}.D does not exist.");
	   if (!mxIsNumeric(BlockD))
	      mexErrMsgTxt("Field BiD{i}.D must be numerical.");
	   pr=(double *)mxGetPr(BlockD);
	   BiDR.valD[i]=pr;
	   if (k==1 && !invert_blocks) {
	      if (spd)
		 pr[0]=1.0/(pr[0]*pr[0]);
	      else
		 pr[0]=1.0/pr[0];
	   }
	 
	   Block=mxGetCell(BL_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BL{i} must be a structure.");

	   BlockI=mxGetField(Block,0,"I");
	   if (BlockI==NULL)
	      mexErrMsgTxt("Field BL{i}.I does not exist.");
	   if (!mxIsNumeric(BlockI))
	      mexErrMsgTxt("Field BL{i}.I must be numerical.");
	   l=mxGetN(BlockI)*mxGetM(BlockI);
	   BLR.nblockrow[i]=l;
	   BLR.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
	   p=BLR.rowind[i];
	   pr=(double *)mxGetPr(BlockI);
	   /* copy and shift indices */
	   for (j=0; j<l; j++)
	       *p++=*pr++-1;
	   
	   BlockE=mxGetField(Block,0,"L");
	   if (BlockE==NULL)
	      mexErrMsgTxt("Field BL{i}.L does not exist.");
	   if (!mxIsNumeric(BlockE))
	      mexErrMsgTxt("Field BL{i}.L must be numerical.");
	   pr=(double *)mxGetPr(BlockE);
	   BLR.valE[i]=pr;
	   	 
       } /* end for i */
#ifdef PRINT_INFO
       mexPrintf("bildlsol: data structures converted, n=%ld,m=%ld\n",n,m);
       fflush(stdout);
#endif

       /* pointer to real input data vector */
       pr=mxGetPr(y_input);
       /* create double real output vector */
       plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)m, mxREAL);
       x_output =plhs[0];
       /* pointer to real output data vector */
       rhssolr=mxGetPr(x_output);
       /* copy rhs to sol */
       if (!mxIsSparse(y_input))
	  memcpy(rhssolr, pr, (size_t)n*m*sizeof(double));
       else {
	  /* init rhssol with 0 */
	  for (i=0; i<n*m; i++)
	      rhssolr[i]=0.0;
	  /* copy y sparse -> full rhs */
	  ja     =(mwIndex *)mxGetIr(y_input);
	  ia     =(mwIndex *)mxGetJc(y_input);
	  for (i=0; i<m; i++) {
	      for (j=ia[i]; j<ia[i+1]; j++) {
		  k=ja[j];
		  rval=pr[j];
		  rhssolr[k+i*n]=rval;
	      } /* end for j */
	  } /* end for i */
       } /* end if */

#ifdef PRINT_INFO
       mexPrintf("bildlsol: call Dbildlsol(%d)\n",++cnt_bildl);
       /*
       for (i=0; i<n; i++)
	 mexPrintf("%8.1le",rhssolr[i]);
       mexPrintf("\n");
       */
       fflush(stdout);
#endif
       /* solve triangular system */
       DSYMbilusol(&BLR,&BiDR,&BLR,pivots, rhssolr,NULL,n,m);
       if (npivots>0)
	  free(pivots);
#ifdef PRINT_INFO
       mexPrintf("bildlsol: Dbildlsol completed\n");
       /*
       for (i=0; i<n; i++)
	 mexPrintf("%8.1le",rhssolr[i]);
       mexPrintf("\n");
       */
       fflush(stdout);
#endif

       for (i=0; i<nblocks; i++) {
	   free(BLR.rowind[i]);
	   
	   Block=mxGetCell(BiD_input,i);
	   BlockJ=mxGetField(Block,0,"J");
	   k=mxGetN(BlockJ)*mxGetM(BlockJ);
	   if (k==1 && !invert_blocks) {
	      BlockD=mxGetField(Block,0,"D");
	      pr=(double *)mxGetPr(BlockD);
	      if (spd)
		 pr[0]=1.0/sqrt(pr[0]);
	      else
		 pr[0]=1.0/pr[0];
	   }
       }
       free(BLR.nblockrow);
       free(BLR.rowind);
       free(BLR.valE);

       free(BiDR.nblockcol);
       free(BiDR.valD);
#ifdef PRINT_INFO
    mexPrintf("bildlsol: memory released\n");fflush(stdout);
#endif
    }
    else {
       /* transform data structures and right hand sides */
       BiDC.nblocks=nblocks;
       BiDC.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
       BiDC.nblockrow=NULL;
       BiDC.colind   =NULL;
       BiDC.rowind   =NULL;
       BiDC.valD     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
       BiDC.valE=NULL;

       BLC.nblocks=nblocks;
       BLC.nblockcol=NULL;
       BLC.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
       BLC.colind   =NULL;
       BLC.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
       BLC.valE     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
       BLC.valD=NULL;

       
       /* map input cell array to struct */
       for (i=0; i<nblocks; i++) {
	   Block=mxGetCell(BiD_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BiD{i} must be a structure.");
	   
	   BlockJ=mxGetField(Block,0,"J");
	   if (BlockJ==NULL)
	      mexErrMsgTxt("Field BiD{i}.J does not exist.");
	   if (!mxIsNumeric(BlockJ))
	      mexErrMsgTxt("Field BiD{i}.J must be numerical.");
	   k=mxGetN(BlockJ)*mxGetM(BlockJ);
	   BiDC.nblockcol[i]=k;

	   BlockD=mxGetField(Block,0,"D");
	   if (BlockD==NULL)
	      mexErrMsgTxt("Field BiD{i}.D does not exist.");
	   if (!mxIsNumeric(BlockD))
	      mexErrMsgTxt("Field BiD{i}.D must be numerical.");
	   BiDC.valD[i]=(doublecomplex *)malloc((size_t)k*k*sizeof(doublecomplex));
	   pD=BiDC.valD[i];
	   
	   pr=(double *)mxGetPr(BlockD);
	   if (mxIsComplex(BlockD)) {
	      pi=(double *)mxGetPi(BlockD);
	      for (j=0; j<k*k; j++,pD++) {
		  pD->r=*pr++;
		  pD->i=*pi++;
	      }
	   }
	   else {
	      for (j=0; j<k*k; j++,pD++) {
		  pD->r=*pr++;
		  pD->i=0.0;
	      }
	   }
	   if (k==1 && !invert_blocks) {
	      pD=BiDC.valD[i];
	      if (spd)
		 pD->r= 1.0/sqrt(pD->r);
	      else {
		 rval=sqrt(pD->r*pD->r+pD->i*pD->i);
		 pD->r= pD->r/rval;
		 pD->i=-pD->i/rval;
	      }
	   }
	   
	 
	   Block=mxGetCell(BL_input,i);
	   if (!mxIsStruct(Block))
	      mexErrMsgTxt("Field BL{i} must be a structure.");

	   BlockI=mxGetField(Block,0,"I");
	   if (BlockI==NULL)
	      mexErrMsgTxt("Field BL{i}.I does not exist.");
	   if (!mxIsNumeric(BlockI))
	      mexErrMsgTxt("Field BL{i}.I must be numerical.");
	   l=mxGetN(BlockI)*mxGetM(BlockI);
	   BLC.nblockrow[i]=l;
	   BLC.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
	   p=BLC.rowind[i];
	   pr=(double *)mxGetPr(BlockI);
	   /* copy and shift indices */
	   for (j=0; j<l; j++)
	       *p++=*pr++-1;
	   
	   BlockE=mxGetField(Block,0,"L");
	   if (BlockE==NULL)
	      mexErrMsgTxt("Field BL{i}.L does not exist.");
	   if (!mxIsNumeric(BlockE))
	      mexErrMsgTxt("Field BL{i}.L must be numerical.");
	   BLC.valE[i]=(doublecomplex *)malloc((size_t)l*k*sizeof(doublecomplex));
	   pE=BLC.valE[i];
	   
	   pr=(double *)mxGetPr(BlockE);
	   if (mxIsComplex(BlockE)) {
	      pi=(double *)mxGetPi(BlockE);
	      for (j=0; j<l*k; j++,pE++) {
		  pE->r=*pr++;
		  pE->i=*pi++;
	      }
	   }
	   else {
	      for (j=0; j<l*k; j++,pE++) {
		  pE->r=*pr++;
		  pE->i=0.0;
	      }
	   }
	   
	 
       } /* end for i */
#ifdef PRINT_INFO
       mexPrintf("bildlsol: data structures converted\n");
       fflush(stdout);
#endif

       /* pointer to non-real output data vector */
       rhssolc=(doublecomplex *)malloc((size_t)n*m*sizeof(doublecomplex));
       pz=rhssolc;
       /* copy rhs to sol */
       /* pointer to real input data vector */
       pr=mxGetPr(y_input);
       if (!mxIsSparse(y_input)) {
	  if (mxIsComplex(y_input)) {
	     pi=mxGetPi(y_input);
	     for (i=0; i<n*m; i++,pz++) {
	         pz->r=*pr++;
		 pz->i=*pi++;
	     }
	  }
	  else {
	     for (i=0; i<n*m; i++,pz++) {
	         pz->r=*pr++;
		 pz->i=0.0;
	     }
	  }
       }
       else {
	  /* init rhssol with 0 */
	  for (i=0; i<n*m; i++)
	      rhssolc[i].r=rhssolc[i].i=0.0;
	  /* copy y sparse -> full rhs */
	  ja     =(mwIndex *)mxGetIr(y_input);
	  ia     =(mwIndex *)mxGetJc(y_input);
	  if (mxIsComplex(y_input)) {
	     pi=mxGetPi(y_input);
	     for (i=0; i<m; i++) {
	         for (j=ia[i]; j<ia[i+1]; j++) {
		     k=ja[j];
		     val.r=pr[j];
		     val.i=pi[j];
		     rhssolc[k+i*n]=val;
		 } /* end for j */
	     } /* end for i */
	  } /* end if */
	  else {
	     for (i=0; i<m; i++) {
	         for (j=ia[i]; j<ia[i+1]; j++) {
		     k=ja[j];
		     val.r=pr[j];
		     val.i=0.0;
		     rhssolc[k+i*n]=val;
		 } /* end for j */
	     } /* end for i */
	  } /* end if-else */
       }
       
#ifdef PRINT_INFO
       mexPrintf("bildlsol: call Zbildlsol(%d)\n",++cnt_bildl);
       /*
       for (i=0; i<n; i++)
	   mexPrintf("%8.1le",rhssolc[i].r);
       mexPrintf("\n");
       for (i=0; i<n; i++)
	   mexPrintf("%8.1le",rhssolc[i].i);
       mexPrintf("\n\n");
       */
       fflush(stdout);
#endif
       /* solve triangular system */
       ZHERbilusol(&BLC,&BiDC,&BLC,pivots, rhssolc,NULL,n,m);
       if (npivots>0)
	  free(pivots);
#ifdef PRINT_INFO
       mexPrintf("bildlsol: Zbildlsol completed\n");
       fflush(stdout);
#endif
       /* create double complex output vector */
       plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)m, mxCOMPLEX);
       x_output =plhs[0];
       pr=mxGetPr(x_output);
       pi=mxGetPi(x_output);
       pz=rhssolc;
       for (i=0; i<n*m; i++,pz++) {
	   *pr++=pz->r;
	   *pi++=pz->i;
       }

       /* release complex-specific memory */
       free(BiDC.nblockcol);
       free(BLC.nblockrow);
       for (i=0; i<nblocks; i++) {
	   free(BiDC.valD[i]);
	   free(BLC.rowind[i]);
	   free(BLC.valE[i]);
       }
       free(BiDC.valD);
       free(BLC.rowind);
       free(BLC.valE);

       free(rhssolc);

#ifdef PRINT_INFO
    mexPrintf("bildlsol: memory released\n");fflush(stdout);
#endif
    }
    
    return;       

}
