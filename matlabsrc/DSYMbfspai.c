/* ========================================================================== */
/* === DSYMbfspai mexFunction ============================================== */
/* ========================================================================== */

/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    $Id: DSYMbfspai.c 6259 2020-05-13 18:24:02Z bolle $

    Usage:

    Return selected inverse 'iA'
    
    Example:

    % for initializing parameters
    iA=DSYMbfspai(PREC,tol)

    given an approximate block factorization P^T S A S P ~ L D^{-1} L^T  up to
    some accuracy, compute approximately iA ~ A^{-1} 
    P,S,L,D^{-1} are stored as part of PREC



    Notice:

        JANUS block ILU

	Copyright (c) 2020 by TU Braunschweig.  All Rights Reserved.

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
#include <omp.h>

#include <blas.h>
#include <janus.h>

#define _DOUBLE_REAL_
#include <ilupackmacros.h>



#define MAX(A,B) (((A)>=(B))?(A):(B))
#define MIN(A,B) (((A)>=(B))?(B):(A))
/* #define PRINT_CHECK     */
/* #define PRINT_INFO      */
/* #define PRINT_INFOnnz   */
/* #define PRINT_INFO1     */



/* ========================================================================== */
/* === mexFunction ========================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right-hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
   /* const char     **fnames; */
   mxArray        *BL, *BiD, *P, *SL,
                  *BL_block, *BL_blockI, *BL_blockJ, *BL_blockL,
                  *BiD_block, *BiD_blockD,
                  *PREC_input, *tol_input, *opts_input, *tmp;
   integer        i,j,k,l,n,*ja,ii,jj,*ia,nnz,
                  nblocks, *idxl;
   size_t         mrows, ncols;
   double         *S_valuesR, mytol, *prBLL, *prBiD, *prBLI, *a,
                  *prBLJ,
                  *Ainv_valuesR;
   mwIndex        *P_ja,          /* row indices of input matrix P         */
                  *P_ia,          /* column pointers of input matrix P     */
                  *S_ja,          /* row indices of input matrix S         */
                  *S_ia,          /* column pointers of input matrix S     */
                  *Ainv_ja,       /* row indices of input matrix Ainv      */
                  *Ainv_ia;       /* column pointers of input matrix Ainv  */

   
    if (nrhs!=2)
       mexErrMsgTxt("2 input arguments required.");
    else if (nlhs!=1)
       mexErrMsgTxt("wrong number of output arguments.");
    else if (!mxIsStruct(prhs[0]))
       mexErrMsgTxt("First input must be a structure.");
    else if (!mxIsNumeric(prhs[1]))
       mexErrMsgTxt("Second input must be a number.");     
    
    /* The first input must be a structure containing L matrix.*/
    PREC_input=(mxArray *)prhs[0];
   
    /* The second input must be a number */
    tol_input=(mxArray *)prhs[1];
    /* get size of input tol */
    mrows=mxGetM(tol_input);
    ncols=mxGetN(tol_input);
    if (mrows*ncols!=1) {
       mexErrMsgTxt("Fifth input must be a number.");
    }
    /* mytol=0.1*(*mxGetPr(tol_input)); */
    mytol=*mxGetPr(tol_input);
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter tol imported\n");fflush(stdout);
#endif
    



    /* extract elements from PREC */
    DSparseBlockMatrix DBL,DBiD;
    integer *p, *invq;
    DSparseMatrix iA;


    
    /* 0. n */
    tmp=mxGetField(PREC_input,0,"n");
    if (!mxIsNumeric(tmp))
       mexErrMsgTxt ("PREC.n must be a number!") ;
    n=*mxGetPr(tmp);
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter PREC.n imported\n");fflush(stdout);
#endif
    
    /* 1. BL */
    BL=mxGetField(PREC_input,0,"BL");
    /* get size of input cell array BL */
    if (!mxIsCell(BL))
       mexErrMsgTxt ("PREC.BL must be a cell array!") ;
    nblocks=MAX(mxGetM(BL),mxGetN(BL));
    if (mxGetM(BL)!=1 && mxGetN(BL)!=1)
       mexErrMsgTxt("PREC.BL must be a 1-dim. cell array!");
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter PREC.BL imported\n");fflush(stdout);
#endif
    DBL.nr=DBL.nc=n;
    DBL.nblocks=nblocks;
    DBL.isreal=1;
    DBL.issingle=0;
    
    DBL.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    DBL.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
    DBL.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    DBL.colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    DBL.valD     =NULL;
    DBL.valE     =(double **) malloc((size_t)nblocks*sizeof(double *));
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: PREC.BL prepared\n");fflush(stdout);
#endif
    for (j=0; j<nblocks; j++) {
        BL_block=mxGetCell(BL,j);

	/* copy PREC.BL{j}.J */
	BL_blockJ=mxGetField(BL_block,0,"J");
	ncols=mxGetN(BL_blockJ)*mxGetM(BL_blockJ);
	DBL.nblockcol[j]=ncols;
	DBL.colind[j]=(integer *)malloc((size_t)ncols*sizeof(integer));
	ja=DBL.colind[j];
	prBLJ=(double *)mxGetPr(BL_blockJ);
	for (k=0; k<ncols; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLJ++;
	    ja[k]=i-1;
	} /* end for k */
	
	/* copy PREC.BL{j}.I */
	BL_blockI=mxGetField(BL_block,0,"I");
	mrows=mxGetN(BL_blockI)*mxGetM(BL_blockI);
	DBL.nblockrow[j]=mrows;
	DBL.rowind[j]=(integer *)malloc((size_t)mrows*sizeof(integer));
	ia=DBL.rowind[j];
	prBLI=(double *)mxGetPr(BL_blockI);
	for (k=0; k<mrows; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLI++;
	    ia[k]=i-1;
	} /* end for k */
	
	/* copy PREC.BL{j}.L */
	BL_blockL=mxGetField(BL_block,0,"L");
	prBLL=(double *)mxGetPr(BL_blockL);
	DBL.valE[j]=(double *) malloc((size_t)mrows*ncols*sizeof(double));
	memcpy(DBL.valE[j],prBLL,mrows*ncols*sizeof(double));
    }
    
    /* 2. BiD */
    BiD=mxGetField(PREC_input,0,"BiD");
    if (!mxIsCell(BiD))
       mexErrMsgTxt ("PREC.BiD must be a cell array!") ;
    /* get size of input matrix BiD */
    if (MAX(mxGetM(BiD),mxGetN(BiD))!=nblocks)
       mexErrMsgTxt("PREC.BiD must be a cell array of same size as the first input!");
    if (mxGetM(BiD)!=1 && mxGetN(BiD)!=1)
       mexErrMsgTxt("PREC.BiD must be a 1-dim. cell array!");
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter PREC.BiD imported\n");fflush(stdout);
#endif
    DBiD.nr=DBiD.nc=n;
    DBiD.nblocks=nblocks;
    DBiD.isreal=1;
    DBiD.issingle=0;
    
    DBiD.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    DBiD.nblockrow=NULL;
    DBiD.rowind   =NULL;
    DBiD.colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    DBiD.valD     =(double **) malloc((size_t)nblocks*sizeof(double *));
    DBiD.valE     =NULL;
    for (j=0; j<nblocks; j++) {
        BiD_block=mxGetCell(BiD,j);

	/* copy PREC.BiD{j}.J */
	BL_blockJ=mxGetField(BiD_block,0,"J");
	ncols=mxGetN(BL_blockJ)*mxGetM(BL_blockJ);
	DBiD.nblockcol[j]=ncols;
	DBiD.colind[j]=(integer *)malloc((size_t)ncols*sizeof(integer));
	ja=DBiD.colind[j];
	prBLJ=(double *)mxGetPr(BL_blockJ);
	for (k=0; k<ncols; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLJ++;
	    ja[k]=i-1;
	} /* end for k */
		
	/* copy PREC.BiD{j}.D */
	BiD_blockD=mxGetField(BiD_block,0,"D");
	prBiD=(double *)mxGetPr(BiD_blockD);
	DBiD.valD[j]  =(double *) malloc((size_t)ncols*ncols*sizeof(double));
	memcpy(DBiD.valD[j],prBiD,ncols*ncols*sizeof(double));
    }

    
    
    /* 3. P */
    P=mxGetField(PREC_input,0,"P");
    if (!mxIsSparse(P))
       mexErrMsgTxt("PREC.P must be in sparse format!") ;
    /* get size of input matrix P */
    mrows=mxGetM(P);
    ncols=mxGetN(P);
    if (mrows!=n || n!=ncols)
       mexErrMsgTxt("PREC.P must be a square matrix of same size as PREC!");
    P_ja=(mwIndex *)mxGetIr(P);
    P_ia=(mwIndex *)mxGetJc(P);
    for (j=0; j<n; j++) {
        if (P_ia[j+1]-P_ia[j]!=1)
	   mexErrMsgTxt("PREC.P must be a permutation matrix!");
    } /* end for j */   
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter PREC.P imported\n");fflush(stdout);
#endif
    p   =(integer *) malloc((size_t)n*sizeof(integer));
    invq=(integer *) malloc((size_t)n*sizeof(integer));
    for (j=0; j<n; j++)
        p[j]=P_ja[j];
    for (j=0; j<n; j++)
        invq[p[j]]=j;
    
    /* 4. SL */
    SL=mxGetField(PREC_input,0,"SL");
    if (!mxIsSparse (SL))
       mexErrMsgTxt ("PREC.SL input matrix must be in sparse format!") ;
    mrows=mxGetM(SL);
    ncols=mxGetN(SL);
    if (mrows!=n || n!=ncols)
       mexErrMsgTxt("PREC.SL must be a square matrix of same size as PREC!");
    S_ja     =(mwIndex *)mxGetIr(SL);
    S_ia     =(mwIndex *)mxGetJc(SL);
    S_valuesR=(double *) mxGetPr(SL);
    for (j=0; j<n; j++) {
        if (S_ia[j+1]-S_ia[j]!=1)
	   mexErrMsgTxt("PREC.SL must be a diagonal matrix!");
	else if (S_ja[j]!=j)
	   mexErrMsgTxt("PREC.SL must be a diagonal matrix!");
    } /* end for j */
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: input parameter PREC.SL imported\n");fflush(stdout);
#endif

    DSYMbfspai(&iA,&DBL,&DBiD,NULL,
	       S_valuesR,NULL,
	       p,invq,mytol);
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: symbfspai completed\n");fflush(stdout);
#endif

    
    /* 4. complete iA with its transpose for the MATLAB interface data */

        
    /* export data structures */
    /* export Ainv */
    nnz=0;
    for (j=0; j<n; j++) 
        nnz+=iA.ncol[j];
    plhs[0]=mxCreateSparse((mwSize)n,(mwSize)n, (mwSize)nnz, mxREAL);
    Ainv_ja         = (mwIndex *)mxGetIr(plhs[0]);
    Ainv_ia         = (mwIndex *)mxGetJc(plhs[0]);
    Ainv_valuesR    = (double *) mxGetPr(plhs[0]);
    Ainv_ia[0]=0;
    for (j=0; j<n; j++) {
        l=iA.ncol[j];
	/* short cuts */
	ja=iA.rowind[j];
	a =iA.val[j];
	k=Ainv_ia[j];
        for (i=0; i<l; i++) {
	    Ainv_ja[k+i]     =ja[i];
	    Ainv_valuesR[k+i]=a[i];
	} /* end for i */
        Ainv_ia[j+1]=Ainv_ia[j]+l;
    } /* end for j */
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: iA exported\n");fflush(stdout);
#endif
        

    /* release memory */
    DSparseDelete(&iA);
    DSparseBlockDelete(&DBL);
    DSparseBlockDelete(&DBiD);
    free(p);
    free(invq);
    
#ifdef PRINT_INFO
    mexPrintf("DSYMbfspai: memory released\n");fflush(stdout);
#endif
    
} /* end DSYMbfspai */


