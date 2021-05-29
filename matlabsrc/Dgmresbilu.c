/* $Id: Dgmresbilu.c 4177 2018-04-30 05:55:52Z bolle $ */
/* ========================================================================== */
/* ==== Dgmresbilu mexFunction ============================================= */
/* ========================================================================== */

/*
    Usage:

    Return computed solution by GMRES with BILU preconditioning
    
    Example:

    % STANDARD GMRES call
    [sol,flag,iter,resvec]=Dgmresbilu(A,rhs,restart,restol,maxit,BL,BiD,BUT,P,Q,SL,SR,pivots,x0);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	April  9, 2017. JANUS Block ILU R1.0

    Notice:

	Copyright (c) 2017 by TU Braunschweig. 
        All Rights Reserved.

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

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

#define MAX(A,B) (((A)>=(B))?(A):(B))
#define MIN(A,B) (((A)>=(B))?(B):(A))

/* #define PRINT_INFO */

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

    const char **fnames;       /* pointers to field names */
    const mwSize  *dims;
    mxClassID  *classIDflags;
    mwSize     mrows, ncols, ndim, mydims[2];
    mwIndex    *A_ia,*A_ja;
    mxArray    *A_input, *rhs_input, *restol_input, *maxit_input, *x0_input,
               *Block, *BlockJ, *BlockI, *BlockE, *BlockD, 
               *restart_input, *BUT_input, *q_input, *SR_input, 
               *BL_input, *BiD_input,
               *p_input, *SL_input;
    int        ifield, nfields;
    integer    iter, *p,*q, *pI,*invq,nblocks,i,j,k,l,m,n, ipar[16],
               maxit,flag, restart, ierr,*pivots,npivots, invert_blocks;
    double     *A_valuesR, restol, *pr, *SL, *SR, nrm1, nrminf, locnrm,
               fpar[16], *w, *y, *sol, *resvec, *rhs,v;
#ifndef _USE_MKL_
    DSparseMatrix A;
#endif
    DSparseBlockMatrix BL,BiD,BUT;

    
    if (nrhs!=14)
       mexErrMsgTxt("Dgmresbilu: fourteen input arguments required.");
    else if (nlhs!=4)
       mexErrMsgTxt("Dgmresbilu: four output arguments are required.");



    /* The first input must be a sparse square matrix */
    A_input=(mxArray *)prhs[0];
    if (!mxIsSparse(A_input))
       mexErrMsgTxt("Dgmresbilu: first input must be a sparse matrix.");
    mrows=mxGetM(A_input);
    ncols=mxGetN(A_input);
    if (mrows!=ncols) {
       mexErrMsgTxt("Dgmresbilu: first input must be a square matrix.");
    }
    A_ja     =(mwIndex *)mxGetIr(A_input);
    A_ia     =(mwIndex *)mxGetJc(A_input);
    A_valuesR=(double *) mxGetPr(A_input);
    n=ncols;
#ifndef _USE_MKL_
    /* convert MATLAB's format to SparseMatrix format */
    A.nr=A.nc=n;
    A.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
    A.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
    A.val   =(double **) malloc((size_t)n*sizeof(double *));
    for (i=0; i<n; i++) {
        A.ncol[i]  =A_ia[i+1]-A_ia[i];
        A.rowind[i]=A_ja     +A_ia[i];
	A.val[i]   =A_valuesR+A_ia[i];
    } /* end for i */
#endif
#ifdef PRINT_INFO
    mexPrintf("sparse matrix A imported\n");fflush(stdout);
#endif
    
    /* Get second input argument */
    rhs_input=(mxArray *)prhs[1];
    if (!mxIsNumeric(rhs_input))
       mexErrMsgTxt("Dgmresbilu: second input must be a vector.");
    mrows=mxGetM(rhs_input);
    ncols=mxGetN(rhs_input);
    if (mrows!=n || ncols!=1) {
       mexErrMsgTxt("Dgmresbilu: second input must be a vector.");
    }
    rhs=mxGetPr(rhs_input);
#ifdef PRINT_INFO
    mexPrintf("right hand side imported\n");fflush(stdout);
#endif
    
    /* Get third input argument */
    restart_input=(mxArray *)prhs[2];
    if (!mxIsNumeric(restart_input))
       mexErrMsgTxt("Dgmresbilu: third input must be a numeric.");
    mrows=mxGetM(restart_input);
    ncols=mxGetN(restart_input);
    if (mrows!=1 || ncols!=1) {
       mexErrMsgTxt("Dgmresbilu: third input must be a scalar.");
    }
    restart=*mxGetPr(restart_input);
#ifdef PRINT_INFO
    mexPrintf("restart imported\n");fflush(stdout);
#endif
    
    /* Get fourth input argument */
    restol_input=(mxArray *)prhs[3];
    if (!mxIsNumeric(restol_input))
       mexErrMsgTxt("Dgmresbilu: fourth input must be a numeric.");
    mrows=mxGetM(restol_input);
    ncols=mxGetN(restol_input);
    if (mrows!=1 || ncols!=1) {
       mexErrMsgTxt("Dgmresbilu: fourth input must be a scalar.");
    }
    restol=*mxGetPr(restol_input);
#ifdef PRINT_INFO
    mexPrintf("restol imported\n");fflush(stdout);
#endif
    
    /* Get fifth input argument */
    maxit_input=(mxArray *)prhs[4];
    if (!mxIsNumeric(maxit_input))
       mexErrMsgTxt("Dgmresbilu: fifth input must be a numeric");
    mrows=mxGetM(maxit_input);
    ncols=mxGetN(maxit_input);
    if (mrows!=1 || ncols!=1) {
       mexErrMsgTxt("Dgmresbilu: fifth input must be a scalar.");
    }
    maxit=*mxGetPr(maxit_input);
#ifdef PRINT_INFO
    mexPrintf("maxit imported\n");fflush(stdout);
#endif
    
    /* The sixth input must be a cell array */
    BL_input=(mxArray *)prhs[5];
    if (!mxIsCell(BL_input)) {
       mexErrMsgTxt("Dgmresbilu: sixth input matrix must be in cell format.") ;
    }
    /* get size of input cell array BL */
    nblocks=mxGetM(BL_input)*mxGetN(BL_input);
#ifdef PRINT_INFO
    mexPrintf("BL imported\n");fflush(stdout);
#endif
    
    /* The sevenh input must be a cell array */
    BiD_input=(mxArray *)prhs[6];
    if (!mxIsCell(BiD_input)) {
       mexErrMsgTxt("Dgmresbilu: seventh input matrix must be in cell format.") ;
    }
#ifdef PRINT_INFO
    mexPrintf("BiD imported\n");fflush(stdout);
#endif

    /* The eighth input must be a cell array */
    BUT_input=(mxArray *)prhs[7];
    if (!mxIsCell(BUT_input)) {
       mexErrMsgTxt("Dgmresbilu: eighth input matrix must be in cell format.") ;
    }
#ifdef PRINT_INFO
    mexPrintf("BUT imported\n");fflush(stdout);
#endif
    
    /* The nineth input must be a sparse matrix */
    p_input=(mxArray *)prhs[8];
    if (!mxIsSparse(p_input)) {
       mexErrMsgTxt("Dgmresbilu: nineth input matrix must be a sparse matrix.") ;
    }
    /* extract permutation from p_input */
    p=(mwIndex *)mxGetIr(p_input);
#ifdef PRINT_INFO
    mexPrintf("p imported\n");fflush(stdout);
#endif
    
    /* The tenth input must be a sparse matrix */
    q_input=(mxArray *)prhs[9];
    if (!mxIsSparse(q_input)) {
       mexErrMsgTxt("Dgmresbilu: tenth input matrix must be a sparse matrix.") ;
    }
    /* extract permutation from q_input */
    q=(mwIndex *)mxGetIr(q_input);
    invq=(mwIndex *)malloc((size_t)n*sizeof(mwIndex));
    for (i=0; i<n; i++)
        invq[q[i]]=i;
#ifdef PRINT_INFO
    mexPrintf("q imported\n");fflush(stdout);
#endif
    
    /* The eleventh input must be a sparse matrix */
    SL_input=(mxArray *)prhs[10];
    if (!mxIsSparse(SL_input)) {
       mexErrMsgTxt("Dgmresbilu: eleventh input matrix must be a sparse matrix.") ;
    }
    /* extract diagonal entries from SL_input */
    SL=(double *) mxGetPr(SL_input);
#ifdef PRINT_INFO
    mexPrintf("SL imported\n");fflush(stdout);
#endif
    
    /* The twelveth input must be a sparse matrix */
    SR_input=(mxArray *)prhs[11];
    if (!mxIsSparse(SR_input)) {
       mexErrMsgTxt("Dgmresbilu: twelveth input matrix must be a sparse matrix.") ;
    }
    /* extract diagonal entries from SR_input */
    SR=(double *) mxGetPr(SR_input);
#ifdef PRINT_INFO
    mexPrintf("SR imported\n");fflush(stdout);
#endif
    
    /* The thirteenth input must be pivots */
    if (!mxIsNumeric((mxArray *)prhs[12]))
       mexErrMsgTxt("Dgmresbilu: thirteenth input argument must be numeric.") ;
    npivots=mxGetM((mxArray *)prhs[12])*mxGetN((mxArray *)prhs[12]);
    if (npivots==0) {
       pivots=NULL;
       invert_blocks=1;
    }
    else {
       pivots=(integer *)malloc((size_t)npivots*sizeof(integer));
       pr=mxGetPr((mxArray *)prhs[12]);
       for (i=0; i<npivots; i++)
	   pivots[i]=pr[i];
       invert_blocks=0;
    }
    
    /* The fourteenth input must be a vector */
    x0_input=(mxArray *)prhs[13];
    if (!mxIsNumeric(x0_input))
       mexErrMsgTxt("Dgmresbilu: fourteenth input must be a vector");
    mrows=mxGetM(x0_input);
    ncols=mxGetN(x0_input);
    if (mrows!=n || ncols!=1) {
       mexErrMsgTxt("Dgmresbilu: fourteenth input must be a vector.");
    }
#ifdef PRINT_INFO
    mexPrintf("x0 imported\n");fflush(stdout);
#endif

        
    /* transform data structures and right hand sides */
    BiD.nblocks=nblocks;
    BiD.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    BiD.nblockrow=NULL;
    BiD.colind   =NULL;
    BiD.rowind   =NULL;
    BiD.valD     =(double **) malloc((size_t)nblocks*sizeof(double *));
    BiD.valE=NULL;

    BL.nblocks=nblocks;
    BL.nblockcol=NULL;
    BL.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
    BL.colind   =NULL;
    BL.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    BL.valE     =(double **) malloc((size_t)nblocks*sizeof(double *));
    BL.valD=NULL;

    BUT.nblocks=nblocks;
    BUT.nblockcol=NULL;
    BUT.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
    BUT.colind   =NULL;
    BUT.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    BUT.valE     =(double **) malloc((size_t)nblocks*sizeof(double *));
    BUT.valD=NULL;
       
    /* map input cell array to struct */
    for (i=0; i<nblocks; i++) {
        /* BiD{i} */
        Block=mxGetCell(BiD_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("Dgmresbilu: field BiD{i} must be a structure.");
	   
        /* BiD{i}.J */
	BlockJ=mxGetField(Block,0,"J");
	if (BlockJ==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BiD{i}.J does not exist.");
	if (!mxIsNumeric(BlockJ))
	   mexErrMsgTxt("Dgmresbilu: field BiD{i}.J must be numerical.");
	k=mxGetN(BlockJ)*mxGetM(BlockJ);
	BiD.nblockcol[i]=k;

        /* BiD{i}.D */
	BlockD=mxGetField(Block,0,"D");
	if (BlockD==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BiD{i}.D does not exist.");
	if (!mxIsNumeric(BlockD))
	   mexErrMsgTxt("Dgmresbilu: field BiD{i}.D must be numerical.");
	pr=(double *)mxGetPr(BlockD);
	BiD.valD[i]=pr;
	if (k==1 && !invert_blocks)
	   pr[0]=1.0/pr[0];

	
        /* BL{i} */
	Block=mxGetCell(BL_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("Dgmresbilu: field BL{i} must be a structure.");

        /* BL{i}.I */
	BlockI=mxGetField(Block,0,"I");
	if (BlockI==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BL{i}.I does not exist.");
	if (!mxIsNumeric(BlockI))
	   mexErrMsgTxt("Dgmresbilu: field BL{i}.I must be numerical.");
	l=mxGetN(BlockI)*mxGetM(BlockI);
	BL.nblockrow[i]=l;
	BL.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
	pI=BL.rowind[i];
	pr=(double *)mxGetPr(BlockI);
	/* copy and shift indices */
	for (j=0; j<l; j++)
	    *pI++=*pr++-1;
	
        /* BL{i}.L */
	BlockE=mxGetField(Block,0,"L");
	if (BlockE==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BL{i}.L does not exist.");
	if (!mxIsNumeric(BlockE))
	   mexErrMsgTxt("Dgmresbilu: field BL{i}.L must be numerical.");
	pr=(double *)mxGetPr(BlockE);
	BL.valE[i]=pr;

	
        /* BUT{i} */
	Block=mxGetCell(BUT_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("Dgmresbilu: field BUT{i} must be a structure.");

        /* BUT{i}.I */
	BlockI=mxGetField(Block,0,"I");
	if (BlockI==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BUT{i}.I does not exist.");
	if (!mxIsNumeric(BlockI))
	   mexErrMsgTxt("Dgmresbilu: field BUT{i}.I must be numerical.");
	l=mxGetN(BlockI)*mxGetM(BlockI);
	BUT.nblockrow[i]=l;
	BUT.rowind[i]=(integer *)malloc((size_t)l*sizeof(integer));
	pI=BUT.rowind[i];
	pr=(double *)mxGetPr(BlockI);
	/* copy and shift indices */
	for (j=0; j<l; j++)
	    *pI++=*pr++-1;
	
        /* BUT{i}.L */
	BlockE=mxGetField(Block,0,"L");
	if (BlockE==NULL)
	   mexErrMsgTxt("Dgmresbilu: field BUT{i}.L does not exist.");
	if (!mxIsNumeric(BlockE))
	   mexErrMsgTxt("Dgmresbilu: field BUT{i}.L must be numerical.");
	pr=(double *)mxGetPr(BlockE);
	BUT.valE[i]=pr;	
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("BL, BiD, BUT converted\n");fflush(stdout);
#endif


    /* ------------------------------------------------------------------------
       ----- use BILU within preconditioned GMRES (adapted from SPARSKIT) ----- 
    */
    
    /* GMRES set up
       since I prefer the backward error as stopping criterion, I provide
       an estimate for ||A||_2
       use geometric mean of ||A||_inf and ||A||_1 as estimate
    */
    /* double buffer for work space */
    w   =(double *) malloc((size_t)((n+3)*(restart+2)+((restart+1)*restart)/2)*sizeof(double));
    /* use w as buffer for the matrix inf norm */
    for (i=0; i<n; i++)
        w[i]=0.0;
    nrm1=0.0;
    for (j=0; j<n; j++) {
        locnrm=0.0;
        for (k=A_ia[j]; k<A_ia[j+1]; k++) {
	    /* row index i of A(i,j) */
	    i=A_ja[k];
	    /* numerical value |A(i,j)| */
	    v=fabs(A_valuesR[k]);
 	    /* compute 1-norm of A(:,j) */
	    locnrm+=v;
	    /* contribute to 1-norm of A(i,:) */
	    w[i]+=v;
	} /* end for k */
	/* build matrix 1-norm */
        nrm1=MAX(nrm1,locnrm);
    } /* end for i */
    nrminf=0.0;
    for (i=0; i<n; i++) {
        /* build matrix inf-norm */
        nrminf=MAX(nrminf,w[i]);
    } /* end for i */
    /* end computation of ||A||_1, ||A||_inf */
#ifdef PRINT_INFO
    mexPrintf("1-norm/inf-norm computed\n");fflush(stdout);
#endif

    /* double buffer for preconditioning */
    y   =(double *) malloc((size_t)n*sizeof(double));
    resvec=(double *)malloc((size_t)(maxit+1)*sizeof(double));
    
    /* four output arguments */
    nlhs=4;

    /* solution vector, passed back to MATLAB as first output parameter */ 
    mydims[0]=n;
    mydims[1]=1;
    plhs[0]=mxCreateNumericArray(2, mydims, mxDOUBLE_CLASS, mxREAL);
    sol=(double *)mxGetPr(plhs[0]);
    /* initial solution */
    pr=mxGetPr(x0_input);
    memcpy(sol,pr, n);

    ipar[0]=0;       /* first call of gmres */
    ipar[1]=2;       /* use right preconditioning */
    ipar[2]=3;       /* backward error as stopping criterion */
    ipar[3]=(n+3)*(restart+2)+((restart+1)*restart)/2; /* size of w */
    ipar[4]=restart; /* restart length */
    ipar[5]=maxit;   /* maximum number of iteration steps */
    ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

    fpar[0]=restol;  /* relative error tolerance */
    fpar[11]=sqrt(nrm1*nrminf);/* estimate for the matrix-2 norm */
    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

    /* main loop using reverse communication */
#ifdef PRINT_INFO
    mexPrintf("start GMRES\n");fflush(stdout);
#endif
    flag=-1;
    while (flag) {
          Dgmres(&n,rhs,sol,ipar,fpar,w);

	  /* what is required by the solver? */
          w--;
	  /* mexPrintf("ipar[0]=%ld, ipar[6]=%ld, fpar[3]=%8.1le, fpar[4]=%8.1le\n",ipar[0],ipar[6],fpar[3],fpar[4]); */
          switch (ipar[0]) {
          case  1:   /* apply matrix vector multiplication */
          case  2:
#ifdef PRINT_INFO
	        mexPrintf("call for mat-vec\n");fflush(stdout);
#endif
#ifdef _USE_MKL_
	        /* y := A^T*x, since MKL supports CSR instead of CSC
		   void mkl_cspblas_dcsrgemv (const char *transa , 
		                              const MKL_INT *m , 
					      const double *a , 
					      const MKL_INT *ia , 
					      const MKL_INT *ja , 
					      const double *x , 
					      double *y );    
		*/
	        mkl_cspblas_dcsrgemv('t', &n,A_valuesR,A_ia,A_ja,
				     w+ipar[7],w+ipar[8]);
#else
	        DMatVec(&A,w+ipar[7],w+ipar[8]);
#endif
		iter=ipar[6];
		if (iter>0)
		   resvec[iter-1]=fpar[4];
		break;
	  case  3:   /* apply preconditioner */
	  case  4:   
	  case  5:   
	  case  6:   
#ifdef PRINT_INFO
	        mexPrintf("call BILU prec.\n");fflush(stdout);
#endif
	        /* left scaling + permutation */
	        for (i=0; i<n; i++) {
		    k=invq[i];  
		    y[k]=SL[i]*w[i+ipar[7]];
		} /* end for i */
	        Dbilusol(&BL,&BiD,&BUT,pivots, y,w+ipar[8], n,1);
	        /* inverse permutation + right scaling */
		for (i=0; i<n; i++) {
		    k=p[i];  
		    w[k+ipar[8]]=y[i];
		    w[k+ipar[8]]*=SR[k];
		} /* end for i */
		break;
          default:   /* the iterative solver terminated with code=ipar[0] */
                ierr=ipar[0];
                flag=0;
                break;
          } /* end switch */
          w++;
    } /* end while */
#ifdef PRINT_INFO
    mexPrintf("GMRES completed\n");fflush(stdout);
#endif

    
    if (ierr<=1)
       flag=-ierr;
    else if (ierr==5)
       flag=-2;
    else if (ierr==-1)
       flag=1;
    else if (ierr==-2)
       flag=2;
    else if (ierr>0)
       flag=-ierr;
    else 
       mexErrMsgTxt("Dgmresbilu: undefined GMRES error code");

    
    /* why did GMRES stop? */
    /* error? */
    if (ierr) {
       if (ierr==-1)
	  mexPrintf("Dgmresbilu: too many iteration steps\n");
       else if (ierr==-2)
	  mexPrintf("Dgmresbilu: not enough work space\n");
       else if (ierr==-3)
	  mexPrintf("Dgmresbilu: algorithm breaks down\n");
       else 
	  mexPrintf("Dgmresbilu: algorithm stops with error code %ld\n",ierr);
    }
    else { /* success */
       if (ipar[6]>=maxit)
	  mexPrintf("Dgmresbilu: preconditioned GMRES stopped after %ld iteration steps\n", ipar[6]);
       /*
       else
	  mexPrintf("Dgmresbilu: preconditioned GMRES successfully completed with %ld iteration steps\n", ipar[6]);
       */
    }
    /* ------------------- END GMRES preconditioned by BILU ------------------
       ----------------------------------------------------------------------- */

    /* second left hand side: flag */
    plhs[1] =mxCreateDoubleMatrix(1,1, mxREAL);
    pr = (double *) mxGetPr(plhs[1]);
    *pr=flag;

    /* third left hand side: iter */
    plhs[2] =mxCreateDoubleMatrix(1,1, mxREAL);
    pr = (double *) mxGetPr(plhs[2]);
    *pr=iter;
    
    /* fourth left hand side: resvec */
    plhs[3] =mxCreateDoubleMatrix(iter,1, mxREAL);
    pr = (double *) mxGetPr(plhs[3]);
    for (i=0; i<iter; i++)
        pr[i]=resvec[i];
#ifdef PRINT_INFO
    mexPrintf("parameters exported\n");fflush(stdout);
#endif

    
    
    if (npivots>0)
       free(pivots);
    free(w);
    free(y);
    free(invq);
    free(resvec);
#ifndef _USE_MKL_
    free(A.ncol);
    free(A.rowind);
    free(A.val);
#endif
    for (i=0; i<nblocks; i++) {
        free(BL.rowind[i]);
        free(BUT.rowind[i]);

	Block=mxGetCell(BiD_input,i);
	BlockJ=mxGetField(Block,0,"J");
	k=mxGetN(BlockJ)*mxGetM(BlockJ);
	if (k==1 && !invert_blocks) {
	   BlockD=mxGetField(Block,0,"D");
	   pr=(double *)mxGetPr(BlockD);
	   pr[0]=1.0/pr[0];
	}
    }
    free(BiD.nblockcol);
    free(BiD.valD);
    free(BL.nblockrow);
    free(BL.rowind);
    free(BL.valE);
    free(BUT.nblockrow);
    free(BUT.rowind);
    free(BUT.valE);
#ifdef PRINT_INFO
    mexPrintf("memory released\n");fflush(stdout);
#endif



    return;
}

