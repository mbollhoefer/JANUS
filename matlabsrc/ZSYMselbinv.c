/* ========================================================================== */
/* === ZSYMselbinv mexFunction ============================================== */
/* ========================================================================== */

/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    $Id$

    Usage:

    Return selected inverse 'iA'
    
    Example:

    % for initializing parameters
    iA=ZSYMselbinv(PREC)

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
   double         *S_valuesR, mytol, *prBLL,*piBLL, *prBiD,*piBiD,
                  *prBLI, *prBLJ, *Ainv_valuesR,*Ainv_valuesI;
   doublecomplex  *a;
   mwIndex        *P_ja,          /* row indices of input matrix P         */
                  *P_ia,          /* column pointers of input matrix P     */
                  *S_ja,          /* row indices of input matrix S         */
                  *S_ia,          /* column pointers of input matrix S     */
                  *Ainv_ja,       /* row indices of input matrix Ainv      */
                  *Ainv_ia;       /* column pointers of input matrix Ainv  */

   
    if (nrhs!=1)
       mexErrMsgTxt("1 input argument required.");
    else if (nlhs!=1)
       mexErrMsgTxt("wrong number of output arguments.");
    else if (!mxIsStruct(prhs[0]))
       mexErrMsgTxt("First input must be a structure.");
    
    /* The first input must be a structure containing L matrix.*/
    PREC_input=(mxArray *)prhs[0];
   
    



    /* extract elements from PREC */
    ZSparseBlockMatrix ZBL,ZBiD;
    integer *p, *invq;
    ZSparseMatrix iA;


    
    /* 0. n */
    tmp=mxGetField(PREC_input,0,"n");
    if (!mxIsNumeric(tmp))
       mexErrMsgTxt ("PREC.n must be a number!") ;
    n=*mxGetPr(tmp);
#ifdef PRINT_INFO
    mexPrintf("ZSYMselbinv: input parameter PREC.n imported\n");fflush(stdout);
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
    mexPrintf("ZSYMselbinv: input parameter PREC.BL imported\n");fflush(stdout);
#endif
    ZBL.nr=ZBL.nc=n;
    ZBL.nblocks=nblocks;
    ZBL.isreal=0;
    ZBL.issingle=0;
    
    ZBL.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    ZBL.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
    ZBL.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    ZBL.colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    ZBL.valD     =NULL;
    ZBL.valE     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
#ifdef PRINT_INFO
    mexPrintf("ZSYMselbinv: PREC.BL prepared\n");fflush(stdout);
#endif
    for (j=0; j<nblocks; j++) {
        BL_block=mxGetCell(BL,j);

	/* copy PREC.BL{j}.J */
	BL_blockJ=mxGetField(BL_block,0,"J");
	ncols=mxGetN(BL_blockJ)*mxGetM(BL_blockJ);
	ZBL.nblockcol[j]=ncols;
	ZBL.colind[j]=(integer *)malloc((size_t)ncols*sizeof(integer));
	ja=ZBL.colind[j];
	prBLJ=(double *)mxGetPr(BL_blockJ);
	for (k=0; k<ncols; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLJ++;
	    ja[k]=i-1;
	} /* end for k */
	
	/* copy PREC.BL{j}.I */
	BL_blockI=mxGetField(BL_block,0,"I");
	mrows=mxGetN(BL_blockI)*mxGetM(BL_blockI);
	ZBL.nblockrow[j]=mrows;
	ZBL.rowind[j]=(integer *)malloc((size_t)mrows*sizeof(integer));
	ia=ZBL.rowind[j];
	prBLI=(double *)mxGetPr(BL_blockI);
	for (k=0; k<mrows; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLI++;
	    ia[k]=i-1;
	} /* end for k */
	
	/* copy PREC.BL{j}.L */
	BL_blockL=mxGetField(BL_block,0,"L");
	prBLL=(double *)mxGetPr(BL_blockL);
	piBLL=(double *)mxGetPi(BL_blockL);
	ZBL.valE[j]=(doublecomplex *) malloc((size_t)mrows*ncols*sizeof(doublecomplex));
	for (k=0; k<mrows*ncols; k++) 
	    ZBL.valE[j][k].r=*prBLL++;
	if (piBLL!=NULL)
	   for (k=0; k<mrows*ncols; k++)
	       ZBL.valE[j][k].i=*piBLL++;
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
    mexPrintf("ZSYMselbinv: input parameter PREC.BiD imported\n");fflush(stdout);
#endif
    ZBiD.nr=ZBiD.nc=n;
    ZBiD.nblocks=nblocks;
    ZBiD.isreal=0;
    ZBiD.issingle=0;
    
    ZBiD.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    ZBiD.nblockrow=NULL;
    ZBiD.rowind   =NULL;
    ZBiD.colind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    ZBiD.valD     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
    ZBiD.valE     =NULL;
    for (j=0; j<nblocks; j++) {
        BiD_block=mxGetCell(BiD,j);

	/* copy PREC.BiD{j}.J */
	BL_blockJ=mxGetField(BiD_block,0,"J");
	ncols=mxGetN(BL_blockJ)*mxGetM(BL_blockJ);
	ZBiD.nblockcol[j]=ncols;
	ZBiD.colind[j]=(integer *)malloc((size_t)ncols*sizeof(integer));
	ja=ZBiD.colind[j];
	prBLJ=(double *)mxGetPr(BL_blockJ);
	for (k=0; k<ncols; k++) {
	    /* recall that the indices in MATLAB structures start with 1 rather than 0 */
	    i=*prBLJ++;
	    ja[k]=i-1;
	} /* end for k */
		
	/* copy PREC.BiD{j}.D */
	BiD_blockD=mxGetField(BiD_block,0,"D");
	prBiD=(double *)mxGetPr(BiD_blockD);
	piBiD=(double *)mxGetPi(BiD_blockD);
	ZBiD.valD[j]  =(doublecomplex *) malloc((size_t)ncols*ncols*sizeof(doublecomplex));
	for (k=0; k<ncols*ncols; k++) 
	    ZBiD.valD[j][k].r=*prBiD++;
	if (piBiD!=NULL)
	   for (k=0; k<ncols*ncols; k++)
	       ZBiD.valD[j][k].i=*piBiD++;
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
    mexPrintf("ZSYMselbinv: input parameter PREC.P imported\n");fflush(stdout);
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
    mexPrintf("ZSYMselbinv: input parameter PREC.SL imported\n");fflush(stdout);
#endif

    ZSYMselbinv(&iA,&ZBL,&ZBiD,NULL,
		S_valuesR,NULL,
		p,invq);
#ifdef PRINT_INFO
    mexPrintf("ZSYMselbinv: symselbinv completed\n");fflush(stdout);
#endif

    
    /* 4. complete iA with its transpose for the MATLAB interface data */

        
    /* export data structures */
    /* export Ainv */
    nnz=0;
    for (j=0; j<n; j++) 
        nnz+=iA.ncol[j];
    plhs[0]=mxCreateSparse((mwSize)n,(mwSize)n, (mwSize)nnz, mxCOMPLEX);
    Ainv_ja         = (mwIndex *)mxGetIr(plhs[0]);
    Ainv_ia         = (mwIndex *)mxGetJc(plhs[0]);
    Ainv_valuesR    = (double *) mxGetPr(plhs[0]);
    Ainv_valuesI    = (double *) mxGetPi(plhs[0]);
    Ainv_ia[0]=0;
    for (j=0; j<n; j++) {
        l=iA.ncol[j];
	/* short cuts */
	ja=iA.rowind[j];
	a =iA.val[j];
	k=Ainv_ia[j];
        for (i=0; i<l; i++) {
	    Ainv_ja[k+i]     =ja[i];
	    Ainv_valuesR[k+i]=a[i].r;
	    Ainv_valuesI[k+i]=a[i].i;
	} /* end for i */
        Ainv_ia[j+1]=Ainv_ia[j]+l;
    } /* end for j */
#ifdef PRINT_INFO
    mexPrintf("ZSYMselbinv: iA exported\n");fflush(stdout);
#endif
        

    /* release memory */
    ZSparseDelete(&iA);
    ZSparseBlockDelete(&ZBL);
    ZSparseBlockDelete(&ZBiD);
    free(p);
    free(invq);
    
#ifdef PRINT_INFO
    mexPrintf("ZSYMselbinv: memory released\n");fflush(stdout);
#endif
    
} /* end ZSYMselbinv */


