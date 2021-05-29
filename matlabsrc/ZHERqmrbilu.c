/* $Id: ZHERqmrbilu.c 4155 2018-04-24 11:25:39Z bolle $ */
/* ========================================================================== */
/* ==== ZHERqmrbilu mexFunction ============================================= */
/* ========================================================================== */

/*
    Usage:

    Return computed solution by SQMR with BILDL preconditioning
    
    Example:

    % STANDARD simplified QMR call
    [sol,flag,iter,resvec]=ZHERqmrbilu(A,rhs,restol,maxit,BL,BiD,P,SL,pivots,x0);



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

    const char **fnames;       /* pointers to field names */
    const mwSize  *dims;
    mxClassID  *classIDflags;
    mwSize     mrows, ncols, ndim, mydims[2];
    mwIndex    *A_ia,*A_ja;
    mxArray    *A_input, *rhs_input, *restol_input, *maxit_input, *x0_input,
               *Block, *BlockJ, *BlockI, *BlockE, *BlockD, 
               *BL_input, *BiD_input,
               *p_input, *SL_input;
    int        ifield, nfields;
    integer    iter, *p,*pI,*invq,nblocks,i,j,k,l,m,n, ipar[16],maxit,flag, ierr,*pivots,npivots,invert_blocks;
    double     *A_valuesR, *A_valuesI, restol, *pr, *pi, *SL, nrm1inf, locnrm,
               fpar[16], *resvec;
    doublecomplex *zbuff, *rhs, *w, *y, *sol;
    ZSparseMatrix A;
    ZSparseBlockMatrix BL,BiD;

    
    if (nrhs!=10)
       mexErrMsgTxt("ZHERqmrbilu: ten input arguments required.");
    else if (nlhs!=4)
       mexErrMsgTxt("ZHERqmrbilu: four output arguments are required.");



    /* The first input must be a sparse square matrix */
    A_input=(mxArray *)prhs[0];
    if (!mxIsSparse(A_input))
       mexErrMsgTxt("ZHERqmrbilu: first input must be a sparse matrix.");
    mrows=mxGetM(A_input);
    ncols=mxGetN(A_input);
    if (mrows!=ncols) {
       mexErrMsgTxt("ZHERqmrbilu: first input must be a square matrix.");
    }
    A_ja     =(mwIndex *)mxGetIr(A_input);
    A_ia     =(mwIndex *)mxGetJc(A_input);
    A_valuesR=(double *) mxGetPr(A_input);
    A_valuesI=(double *) mxGetPi(A_input);
    n=ncols;
    /* convert MATLAB's format to SparseMatrix format */
    A.nr=A.nc=n;
    A.nnz=A_ia[n];
    A.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
    A.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
    A.val   =(doublecomplex **)malloc((size_t)n*sizeof(doublecomplex *));
    zbuff   =(doublecomplex *) malloc((size_t)A.nnz*sizeof(doublecomplex));
    if (A_valuesI!=NULL)
       for (i=0; i<n; i++) {
	   A.ncol[i]  =A_ia[i+1]-A_ia[i];
	   A.rowind[i]=A_ja     +A_ia[i];
	   A.val[i]   =zbuff    +A_ia[i];
	   for (j=0; j<A.ncol[i]; j++) {
	       A.val[i][j].r=A_valuesR[A_ia[i]+j];
	       A.val[i][j].i=A_valuesI[A_ia[i]+j];
	   } /* end for j */
    } /* end for i */
    else
       for (i=0; i<n; i++) {
	   A.ncol[i]  =A_ia[i+1]-A_ia[i];
	   A.rowind[i]=A_ja     +A_ia[i];
	   A.val[i]   =zbuff    +A_ia[i];
	   for (j=0; j<A.ncol[i]; j++) {
	       A.val[i][j].r=A_valuesR[A_ia[i]+j];
	       A.val[i][j].i=0.0;
	   } /* end for j */
       } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("complex-valued sparse matrix A imported\n");fflush(stdout);
#endif
    
    /* Get second input argument */
    rhs_input=(mxArray *)prhs[1];
    if (!mxIsNumeric(rhs_input))
       mexErrMsgTxt("ZHERqmrbilu: second input must be a vector.");
    mrows=mxGetM(rhs_input);
    ncols=mxGetN(rhs_input);
    if (mrows!=n || ncols!=1) {
       mexErrMsgTxt("ZHERqmrbilu: second input must be a vector.");
    }
    pr=mxGetPr(rhs_input);
    pi=mxGetPi(rhs_input);
    rhs=(doublecomplex *)malloc((size_t)n*sizeof(doublecomplex));
    if (pi!=NULL)
       for (i=0; i<n; i++) {
	   rhs[i].r=*pr++;
	   rhs[i].i=*pi++;
       } /* end for i */
    else
       for (i=0; i<n; i++) {
	   rhs[i].r=*pr++;
	   rhs[i].i=0.0;
       } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("complex-valued right hand side imported\n");fflush(stdout);
#endif
    
    /* Get third input argument */
    restol_input=(mxArray *)prhs[2];
    if (!mxIsNumeric(restol_input))
       mexErrMsgTxt("ZHERqmrbilu: third input must be a numeric.");
    mrows=mxGetM(restol_input);
    ncols=mxGetN(restol_input);
    if (mrows!=1 || ncols!=1) {
       mexErrMsgTxt("ZHERqmrbilu: third input must be a scalar.");
    }
    restol=*mxGetPr(restol_input);
#ifdef PRINT_INFO
    mexPrintf("restol imported\n");fflush(stdout);
#endif
    
    /* Get fourth input argument */
    maxit_input=(mxArray *)prhs[3];
    if (!mxIsNumeric(maxit_input))
       mexErrMsgTxt("ZHERqmrbilu: fourth input must be a numeric");
    mrows=mxGetM(maxit_input);
    ncols=mxGetN(maxit_input);
    if (mrows!=1 || ncols!=1) {
       mexErrMsgTxt("ZHERqmrbilu: fourth input must be a scalar.");
    }
    maxit=*mxGetPr(maxit_input);
#ifdef PRINT_INFO
    mexPrintf("maxit imported\n");fflush(stdout);
#endif
    
    /* The fifth input must be a cell array */
    BL_input=(mxArray *)prhs[4];
    if (!mxIsCell(BL_input)) {
       mexErrMsgTxt("ZHERqmrbilu: fifth input matrix must be in cell format.") ;
    }
    /* get size of input cell array BL */
    nblocks=mxGetM(BL_input)*mxGetN(BL_input);
#ifdef PRINT_INFO
    mexPrintf("BL imported\n");fflush(stdout);
#endif
    
    /* The sixth input must be a cell array */
    BiD_input=(mxArray *)prhs[5];
    if (!mxIsCell(BiD_input)) {
       mexErrMsgTxt("ZHERqmrbilu: sixth input matrix must be in cell format.") ;
    }
#ifdef PRINT_INFO
    mexPrintf("BiD imported\n");fflush(stdout);
#endif

    /* The seventh input must be a sparse matrix */
    p_input=(mxArray *)prhs[6];
    if (!mxIsSparse(p_input)) {
       mexErrMsgTxt("ZHERqmrbilu: seventh input matrix must be a sparse matrix.") ;
    }
    /* extract permutation from p_input */
    p=(mwIndex *)mxGetIr(p_input);
    invq=(mwIndex *)malloc((size_t)n*sizeof(mwIndex));
    for (i=0; i<n; i++)
        invq[p[i]]=i;
#ifdef PRINT_INFO
    mexPrintf("p imported\n");fflush(stdout);
#endif
    
    /* The eighth input must be a sparse matrix */
    SL_input=(mxArray *)prhs[7];
    if (!mxIsSparse(SL_input)) {
       mexErrMsgTxt("ZHERqmrbilu: eighth input matrix must be a sparse matrix.") ;
    }
    /* extract diagonal entries from SL_input */
    SL=(double *) mxGetPr(SL_input);
#ifdef PRINT_INFO
    mexPrintf("SL imported\n");fflush(stdout);
#endif
    
    /* The nineth input must be pivots */
    if (!mxIsNumeric((mxArray *)prhs[8]))
       mexErrMsgTxt("ZHERqmrbilu: nineth input argument must be numeric.") ;
    npivots=mxGetM((mxArray *)prhs[8])*mxGetN((mxArray *)prhs[8]);
    if (npivots==0) {
       pivots=NULL;
       invert_blocks=1;
    }
    else {
       pivots=(integer *)malloc((size_t)n*sizeof(integer));
       pr=mxGetPr((mxArray *)prhs[8]);
       for (i=0; i<npivots; i++)
	   pivots[i]=pr[i];
       invert_blocks=0;
    }
    
    /* The tenth input must be a vector */
    x0_input=(mxArray *)prhs[9];
    if (!mxIsNumeric(x0_input))
       mexErrMsgTxt("ZHERqmrbilu: tenth input must be a vector");
    mrows=mxGetM(x0_input);
    ncols=mxGetN(x0_input);
    if (mrows!=n || ncols!=1) {
       mexErrMsgTxt("ZHERqmrbilu: tenth input must be a vector.");
    }
    /* copy initial solution */
    pr=mxGetPr(x0_input);
    pi=mxGetPi(x0_input);
    sol=(doublecomplex *)malloc((size_t)n*sizeof(doublecomplex));
    if (pi!=NULL)
       for (i=0; i<n; i++) {
	   sol[i].r=*pr++;
	   sol[i].i=*pi++;
       } /* end for i */
    else
       for (i=0; i<n; i++) {
	   sol[i].r=*pr++;
	   sol[i].i=0.0;
       } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("x0 imported\n");fflush(stdout);
#endif


    
    /* transform data structures and right hand sides */
    BiD.nblocks=nblocks;
    BiD.nblockcol=(integer *) malloc((size_t)nblocks*sizeof(integer));
    BiD.nblockrow=NULL;
    BiD.colind   =NULL;
    BiD.rowind   =NULL;
    BiD.valD     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
    BiD.valE=NULL;

    BL.nblocks=nblocks;
    BL.nblockcol=NULL;
    BL.nblockrow=(integer *) malloc((size_t)nblocks*sizeof(integer));
    BL.colind   =NULL;
    BL.rowind   =(integer **)malloc((size_t)nblocks*sizeof(integer *));
    BL.valE     =(doublecomplex **) malloc((size_t)nblocks*sizeof(doublecomplex *));
    BL.valD=NULL;
       
    /* map input cell array to struct */
    for (i=0; i<nblocks; i++) {
        /* BiD{i} */
        Block=mxGetCell(BiD_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("ZHERqmrbilu: field BiD{i} must be a structure.");
	   
        /* BiD{i}.J */
	BlockJ=mxGetField(Block,0,"J");
	if (BlockJ==NULL)
	   mexErrMsgTxt("ZHERqmrbilu: field BiD{i}.J does not exist.");
	if (!mxIsNumeric(BlockJ))
	   mexErrMsgTxt("ZHERqmrbilu: field BiD{i}.J must be numerical.");
	k=mxGetN(BlockJ)*mxGetM(BlockJ);
	BiD.nblockcol[i]=k;

        /* BiD{i}.D */
	BlockD=mxGetField(Block,0,"D");
	if (BlockD==NULL)
	   mexErrMsgTxt("ZHERqmrbilu: field BiD{i}.D does not exist.");
	if (!mxIsNumeric(BlockD))
	   mexErrMsgTxt("ZHERqmrbilu: field BiD{i}.D must be numerical.");
	/* copy data */
	pr=(double *)mxGetPr(BlockD);
	pi=(double *)mxGetPi(BlockD);
	BiD.valD[i]=(doublecomplex *)malloc((size_t)k*k*sizeof(doublecomplex));
	if (pi!=NULL)
	   for (j=0; j<k*k; j++) {
	       BiD.valD[i][j].r=*pr++;
	       BiD.valD[i][j].i=*pi++;
	   } /* end for j */
	else
	   for (j=0; j<k*k; j++) {
	       BiD.valD[i][j].r=*pr++;
	       BiD.valD[i][j].i=0.0;
	   } /* end for j */
	if (k==1 && !invert_blocks) 
	   BiD.valD[i][0].r=1.0/BiD.valD[i][0].r;
	
        /* BL{i} */
	Block=mxGetCell(BL_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("ZHERqmrbilu: field BL{i} must be a structure.");

        /* BL{i}.I */
	BlockI=mxGetField(Block,0,"I");
	if (BlockI==NULL)
	   mexErrMsgTxt("ZHERqmrbilu: field BL{i}.I does not exist.");
	if (!mxIsNumeric(BlockI))
	   mexErrMsgTxt("ZHERqmrbilu: field BL{i}.I must be numerical.");
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
	   mexErrMsgTxt("ZHERqmrbilu: field BL{i}.L does not exist.");
	if (!mxIsNumeric(BlockE))
	   mexErrMsgTxt("ZHERqmrbilu: field BL{i}.L must be numerical.");
	/* copy data */
	pr=(double *)mxGetPr(BlockE);
	pi=(double *)mxGetPi(BlockE);
	BL.valE[i]=(doublecomplex *)malloc((size_t)l*k*sizeof(doublecomplex));
	if (pi!=NULL)
	   for (j=0; j<l*k; j++) {
	       BL.valE[i][j].r=*pr++;
	       BL.valE[i][j].i=*pi++;
	   } /* end for j */
	else
	   for (j=0; j<l*k; j++) {
	       BL.valE[i][j].r=*pr++;
	       BL.valE[i][j].i=0.0;
	   } /* end for j */
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("BL, BiD converted\n");fflush(stdout);
#endif


    /* ------------------------------------------------------------------------
       ----- use BILDL within preconditioned SQMR (adapted from SPARSKIT) ----- 
    */
    
    /* SQMR set up
       since I prefer the backward error as stopping criterion, I provide
       an estimate for ||A||_2 and use ||A||_inf=||A||_1 as estimate */
    
    /* compute 1-norm of A */
    nrm1inf=0.0;
    for (i=0; i<n; i++) {
        locnrm=0.0;
        for (k=A_ia[i]; k<A_ia[i+1]; k++) {
	    /* j=A_ja[k]; */
	    locnrm+=sqrt(A_valuesR[k]*A_valuesR[k]+A_valuesI[k]*A_valuesI[k]);
	} /* end for k */
	nrm1inf=MAX(nrm1inf,locnrm);
    } /* end for i */
    /* end computation of ||A||_1=||A||_inf */
#ifdef PRINT_INFO
    mexPrintf("1-norm computed\n");fflush(stdout);
#endif

    /* double buffer for preconditioning, work space */
    w  =(doublecomplex *)malloc((size_t)7*n*sizeof(doublecomplex));
    y=w+6*n;
    resvec=(double *)malloc((size_t)(maxit+1)*sizeof(double));
    

    ipar[0]=0;       /* first call of sqmr */
    ipar[1]=2;       /* use preconditioning */
    ipar[2]=3;       /* backward error as stopping criterion */
    ipar[3]=n*6;     /* size of w */
    ipar[5]=maxit;   /* maximum number of iteration steps */
    ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

    fpar[0]=restol;  /* relative error tolerance */
    fpar[11]=nrm1inf;/* estimate for the matrix-2 norm */
    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

    /* main loop using reverse communication */
#ifdef PRINT_INFO
    mexPrintf("start SQMR\n");fflush(stdout);
#endif
    flag=-1;
    while (flag) {
          ZHERqmr(&n,rhs,sol,ipar,fpar,w);

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
	        mkl_cspblas_zcsrgemv('t', &n,zbuff,A_ia,A_ja,
				     w+ipar[7],w+ipar[8]);
#else
		/* column-oriented matrix-vector multiplcation */
	        ZMatVec(&A,w+ipar[7],w+ipar[8]);
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
	        mexPrintf("call BILDL prec.\n");fflush(stdout);
#endif
	        /* left scaling + permutation */
	        for (i=0; i<n; i++) {
		    k=invq[i];  
		    y[k].r=SL[i]*w[i+ipar[7]].r;
		    y[k].i=SL[i]*w[i+ipar[7]].i;
		} /* end for i */
	        ZHERbilusol(&BL,&BiD,&BL,pivots, y,w+ipar[8], n,1);
	        /* inverse permutation + right scaling */
		for (i=0; i<n; i++) {
		    k=p[i];  
		    w[k+ipar[8]]=y[i];
		    w[k+ipar[8]].r*=SL[k];
		    w[k+ipar[8]].i*=SL[k];
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
    mexPrintf("SQMR completed\n");fflush(stdout);
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
       mexErrMsgTxt("ZHERqmrbilu: undefined SQMR error code");

    
    /* why did SQMR stop? */
    /* error? */
    if (ierr) {
       if (ierr==-1)
	  mexPrintf("ZHERqmrbilu: too many iteration steps\n");
       else if (ierr==-2)
	  mexPrintf("ZHERqmrbilu: not enough work space\n");
       else if (ierr==-3)
	  mexPrintf("ZHERqmrbilu: algorithm breaks down\n");
       else 
	  mexPrintf("ZHERqmrbilu: algorithm stops with error code %ld\n",ierr);
    }
    else { /* success */
       if (ipar[6]>=maxit)
	  mexPrintf("ZHERqmrbilu: preconditioned SQMR stopped after %ld iteration steps\n", ipar[6]);
       /*
       else
	  mexPrintf("ZHERqmrbilu: preconditioned SQMR successfully completed with %ld iteration steps\n", ipar[6]);
       */
    }
    /* ------------------- END SQMR preconditioned by BILDL ------------------
       ----------------------------------------------------------------------- */

    /* four output arguments */
    nlhs=4;

    /* solution vector, passed back to MATLAB as first output parameter */ 
    mydims[0]=n;
    mydims[1]=1;
    plhs[0]=mxCreateNumericArray(2, mydims, mxDOUBLE_CLASS, mxCOMPLEX);
    /* copy computed solution */
    pr=mxGetPr(plhs[0]);
    pi=mxGetPi(plhs[0]);
    for (i=0; i<n; i++) {
        *pr++=sol[i].r;
	*pi++=sol[i].i;
    } /* end for i */
      
    
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
    free(rhs);
    free(sol);
    free(invq);
    free(resvec);
    free(A.ncol);
    free(A.rowind);
    free(A.val);
    free(zbuff);
    for (i=0; i<nblocks; i++) {
        free(BiD.valD[i]);
        free(BL.rowind[i]);
        free(BL.valE[i]);
    }
    free(BiD.nblockcol);
    free(BiD.valD);
    free(BL.nblockrow);
    free(BL.rowind);
    free(BL.valE);
#ifdef PRINT_INFO
    mexPrintf("memory released\n");fflush(stdout);
#endif



    return;
}

