/* $Id: bilunnz.c 3620 2017-09-09 15:19:03Z bolle $ */
/* ========================================================================== */
/* === bilunnz mexFunction ============================================== */
/* ========================================================================== */

/*
    Usage:

    Given a block factorization BL BiD^{-1} BUT^T, count its fill-in

    Example:

    nz=bilunnz(BL,BiD,BUT)


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	September 22, 2016. JANUS Block ILU R1.0.  

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
    mxArray            *Block, *BlockJ, *BlockI,
                       *BL_input, *BiD_input, *BUT_input, *nz_output;
    double             *pr;
    integer            i,j,k,l,m,nblocks;

    static int cnt_bilu=0;
    
    
    if (nrhs!=3)
       mexErrMsgTxt("Three input arguments required.");
    else if (nlhs!=1)
       mexErrMsgTxt("wrong number of output arguments.");


    /* The first input must be a cell array */
    BL_input=(mxArray *)prhs[0];
    if (!mxIsCell(BL_input)) {
       mexErrMsgTxt("bilunnz: first input matrix must be in cell format.") ;
    }
    /* get size of input cell array BL */
    nblocks=mxGetM(BL_input)*mxGetN(BL_input);

    /* The second input must be a cell array */
    BiD_input=(mxArray *)prhs[1];
    if (!mxIsCell(BiD_input)) {
       mexErrMsgTxt("bilunnz: second input matrix must be in cell format.") ;
    }
    
    /* The third input must be a cell array */
    BUT_input=(mxArray *)prhs[2];
    if (!mxIsCell(BUT_input)) {
       mexErrMsgTxt("bilunnz: third input matrix must be in cell format.") ;
    }


    /* check cell arrays for their sizes */
    j=0;
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

	
	Block=mxGetCell(BL_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("Field BL{i} must be a structure.");

	BlockI=mxGetField(Block,0,"I");
	if (BlockI==NULL)
	   mexErrMsgTxt("Field BL{i}.I does not exist.");
	if (!mxIsNumeric(BlockI))
	   mexErrMsgTxt("Field BL{i}.I must be numerical.");
	l=mxGetN(BlockI)*mxGetM(BlockI);

	
	Block=mxGetCell(BUT_input,i);
	if (!mxIsStruct(Block))
	   mexErrMsgTxt("Field BUT{i} must be a structure.");
	
	BlockI=mxGetField(Block,0,"I");
	if (BlockI==NULL)
	   mexErrMsgTxt("Field BUT{i}.I does not exist.");
	if (!mxIsNumeric(BlockI))
	   mexErrMsgTxt("Field BUT{i}.I must be numerical.");
	m=mxGetN(BlockI)*mxGetM(BlockI);

	j+=k*(k+l+m);
    } /* end for i */
    

    /* create double complex output vector */
    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    nz_output =plhs[0];
    pr=mxGetPr(nz_output);
    *pr=j;
    
    return;       
}
