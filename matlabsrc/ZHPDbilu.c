/* $Id: ZHPDbilu.c 7315 2021-05-28 21:00:20Z bolle $ */
/* ========================================================================== */
/* === ZHPDbilu mexFunction ============================================== */
/* ========================================================================== */

/*
    Usage:

    Given a permutation matrix P as well as a scaling matrix S_L,
    compute a block-structured approximate factorization

          P^T S_L A S_L P ~  BL BiD^{-1} BL^T

    with a block (inverse) diagonal matrix BiD, block unit lower triangular 
    factor BL.


    Example:

    [BL,BiD,P,S_L,determinant,options]=ZHPDbilu(A,options)


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	December 25, 2016. JANUS Block ILU R1.0.  

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

#define MAX_FIELDS 100
#define MAX(A,B) (((A)>=(B))?(A):(B))
#define MIN(A,B) (((A)>=(B))?(B):(A))
#define ELBOW    MAX(4.0,2.0)
#define USE_COSINE_BLOCKS /* */
#define TAU      0.8
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
    const char         *BLnames[]={"J","I", "L","D"},
                       *optionsnames[]={"cosine","blocking_strategy","progressive_aggregation","ordering","droptol","perturbation","invert_blocks","level_of_fill"},
		       *optionsnamesx[]={"cosine","blocking_strategy","progressive_aggregation","ordering","droptol","perturbation","invert_blocks","level_of_fill","blocksize"};
    const char         **fnames;
    const mwSize       *dims;
    char               *pdata, *input_buf, *output_buf;
    mxClassID          *classIDflags;
    mwSize             mydims[1], nnz,ndim, buflen;
    mxArray            *Block, *BlockJ, *BlockI, *BlockL, *BlockD, 
                       *A_input, *options_input,*options_output,
                       *P, *S_L,
                       *BL_output, *BiD_output, 
                       *tmp, *fout;
    mwIndex            *ja,  /* row indices of input matrix A,P,S_L,D */
                       *ia;  /* column pointers of input matrix A,P,S_L,D */
    int                ifield, status, nfields;
    integer            i,j,k,l,ll,jj,m,n,*p,*invq,*invq2,*p2,*invq3,*p3,r,s,mynblocks,*myblocksize,flag=0,
                       nB=0, ierr,bi, *invblock, *blockind, cosine=0, ilu1t=1, pa=1,
                       perturbation=1,
                       invert_blocks=1,
                       level_of_fill=1,
                       cnt, *idx, *idxpos;
    double             *valuesR, *pr, *valuesI, *pi,
                       *pcolscale, *prowscale;
    doublecomplex      *colscale, *rowscale;
    doublecomplex      determinant;
    size_t             mrows, ncols;
    ZILUPACKparam      options;
    Zmat               A,B;
    ZSparseMatrix      C;
    ZSparseBlockMatrix BL,BiD;

    if (nrhs!=2)
       mexErrMsgTxt("Two input arguments required.");
    else if (nlhs!=6)
       mexErrMsgTxt("wrong number of output arguments.");


    /* The first input must be a square matrix.*/
    A_input=(mxArray *)prhs[0];
    /* get size of input matrix A */
    mrows=mxGetM(A_input);
    ncols=mxGetN(A_input);
    nnz  =mxGetNzmax(A_input);
    if (mrows!=ncols) {
       mexErrMsgTxt("First input must be a square matrix.");
    }
    if (!mxIsSparse(A_input)) {
       mexErrMsgTxt("ZHPDbilu: input matrix must be in sparse format.") ;
    }
    n=mrows;

    /* copy input matrix to sparse row format */
    A.nc=A.nr=n;
    A.ia=(integer *)      malloc((size_t)(n+1)*sizeof(integer));
    A.ja=(integer *)      malloc((size_t)nnz  *sizeof(integer));
    A.a =(doublecomplex *)malloc((size_t)nnz  *sizeof(doublecomplex));

    ja     =(mwIndex *)mxGetIr(A_input);
    ia     =(mwIndex *)mxGetJc(A_input);
    valuesR=(double *) mxGetPr(A_input);
    valuesI=(double *) mxGetPi(A_input);

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    
    A.ia[0]=1;
    for (i=0; i<ncols; i++) {
        A.ia[i+1]=A.ia[i];
	for (j=ia[i]; j<ia[i+1]; j++) {
	    k=ja[j];
	    if (k>=i) {
	       l=A.ia[i+1]-1;
	       A.ja[l]=k+1;
	       A.a [l].r=valuesR[j];
	       A.a [l].i=valuesI[j];
	       A.ia[i+1]=l+2;
	    }
	}
    }
    A.nnz=A.ia[n]-1;

    /* initialize ILUPACK options structure to its default options */
    ZHPDAMGinit(&A,&options);
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: sparse matrix imported\n");fflush(stdout);
#endif


    /* Get second input argument */
    options_input=(mxArray *)prhs[1];
    nfields=mxGetNumberOfFields(options_input);

    /* Allocate memory  for storing classIDflags */
    classIDflags=(mxClassID *)mxCalloc((size_t)nfields, (size_t)sizeof(mxClassID));
    
    /* allocate memory  for storing pointers */
    fnames=mxCalloc((size_t)nfields, (size_t)sizeof(*fnames));
    /* Get field name pointers */
    for (ifield=0; ifield<nfields; ifield++) {
        fnames[ifield]=mxGetFieldNameByNumber(options_input,ifield);
    }

    /* import data */
    for (ifield=0; ifield<nfields; ifield++) {
	tmp=mxGetFieldByNumber(options_input,0,ifield);
	classIDflags[ifield]=mxGetClassID(tmp); 

	ndim=mxGetNumberOfDimensions(tmp);
	dims=mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield]==mxCHAR_CLASS) {
	   /* Get the length of the input string. */
	   buflen=(mxGetM(tmp)*mxGetN(tmp))+1;

	   /* Allocate memory for input and output strings. */
	   input_buf=(char *)mxCalloc((size_t)buflen, (size_t)sizeof(char));

	   /* Copy the string data from tmp into a C string input_buf. */
	   status=mxGetString(tmp, input_buf, buflen);
	   
	   /* ordering */
	   if (!strcmp("ordering",fnames[ifield])) {
              if (strcmp(options.ordering,input_buf)) {
		 options.ordering=(char *)malloc((size_t)buflen*sizeof(char));
		 strcpy(options.ordering,input_buf);
		 flag=-1;
	      }
	   }
	   else if (!strcmp("cosine",fnames[ifield])) {
              if (!strcmp("yes",input_buf))
		 cosine=1;
	      else
		 cosine=0;
	   }
           else if (!strcmp("blocking_strategy",fnames[ifield])) {
	      if (!strncasecmp("ilu1t",input_buf,5))
                 ilu1t=BLOCK_ILU1T;
	      else if (!strncasecmp("ilupt",input_buf,5))
	         ilu1t=BLOCK_ILUPT;
	      else if (!strncasecmp("super",input_buf,5))
	         ilu1t=BLOCK_SUPERNODES;
	      else if (!strncasecmp("appsuper",input_buf,8))
	         ilu1t=BLOCK_APP_SUPERNODES;
	      else
	         ilu1t=BLOCK_NONE;
           }
	   else if (!strcmp("progressive_aggregation",fnames[ifield])) {
              if (!strcmp("yes",input_buf))
		 pa=1;
	      else
		 pa=0;
	   }
	   else {
	      /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	   mxFree(input_buf);
	} 
	else {
	   if (!strcmp("droptol",fnames[ifield])) {
	      options.droptol=*mxGetPr(tmp);
	   }
	   else if (!strcmp("perturbation",fnames[ifield])) {
	      if (*mxGetPr(tmp))
		 perturbation=1;
	      else
		 perturbation=0;
	   }
	   else if (!strcmp("invert_blocks",fnames[ifield])) {
	      if (*mxGetPr(tmp))
		 invert_blocks=1;
	      else
		 invert_blocks=0;
	   }
	   else if (!strcmp("level_of_fill",fnames[ifield])) {
	      level_of_fill=*mxGetPr(tmp);
	   }
	   else if (!strcmp("blocksize",fnames[ifield])) {
	      options.nblocks=mxGetM(tmp)*mxGetN(tmp);
	      options.blocksize=(integer *)malloc((size_t)options.nblocks*sizeof(integer));
	      pr=(double *)mxGetPr(tmp);
	      cnt=0;
	      for (i=0; i<options.nblocks; i++) {
		  options.blocksize[i]=pr[i];
		  cnt+=pr[i];
	      } /* end for i */
	      if (cnt!=n) {
		 mexErrMsgTxt("ZHPDbilu: options.blocksize must exactly match the size of the matrix") ;
	      }
	      cnt=0;
	   }
	   else {
	     /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	}
    }
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: options structure imported\n");fflush(stdout);
#endif
    mxFree(classIDflags);
    mxFree(fnames);



    /* memory for initial scaling and permutation */
    p   =(integer *)malloc((size_t)n*sizeof(integer));
    invq=(integer *)malloc((size_t)n*sizeof(integer));
    /* initial permutations in C-style */
    for (i=0; i<n; i++) {
        p[i]=invq[i]=i;
    } /* end for i */
    
    /* switch matrix A back to C notation and adapt data structures for bilu (and cosine_blocks) */
    for (i=0; i<=n; i++)
        A.ia[i]--;
    for (i=0; i<A.ia[n]; i++)
        A.ja[i]--;
    /* use slightly different data structure for bilu and cosine_blocks */
    C.nr=A.nr; C.nc=n;
    C.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
    C.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
    C.val   =(doublecomplex **) malloc((size_t)n*sizeof(doublecomplex *));
    C.nnz=A.nnz;
    for (i=0; i<n; i++) {
        C.ncol[i]=A.ia[i+1]-A.ia[i];
	C.rowind[i]=A.ja+A.ia[i];
	C.val   [i]=A.a +A.ia[i];
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: C set up\n");fflush(stdout);
#endif


    /* if the user does not predefine an initial block structure we do Saad's cosine blocking */
    if (options.blocksize==NULL) {
       myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
       mynblocks=0;
       /* previously computed permutations p and invq will be respected and updated */
#ifdef USE_COSINE_BLOCKS       
       if (cosine) {
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: call modified cosine-based blocking\n");fflush(stdout);
#endif
	  Zmcosine_sblocks(&C, p,invq, myblocksize, &mynblocks, TAU);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: modified cosine-based blocking applied\n");fflush(stdout);
#endif
       }
       else {
	  for (i=0; i<n; i++)
	      myblocksize[i]=1;
	  mynblocks=n;
       } /* end if-else */
#else
       for (i=0; i<n; i++)
	   myblocksize[i]=1;
       mynblocks=n;
#endif
    }
    else {
       myblocksize=options.blocksize;
       mynblocks=options.nblocks;
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: initial permutation turned off\n");fflush(stdout);
#endif
    }
    
    /* row scaling vectors real and complex */
    pcolscale=(double *)malloc((size_t)n*sizeof(double));
    prowscale=pcolscale;
    colscale=(doublecomplex *)malloc((size_t)n*sizeof(doublecomplex));
    rowscale=colscale;

    /* cosine not applied or cosine did not find dense blocks */
    if (mynblocks==n) {
       /* switch from C to FORTRAN style to apply ILUPACK drivers */
       for (i=0; i<=n; i++)
	   A.ia[i]++;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]++;
       /* select permutation and scaling driver by user request */
       if (!strcmp("mtmetis",options.ordering)) {
          ierr=ZHPDperm_mtmetis(A, rowscale, colscale, p, invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: mtmetis performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       else if (!strcmp("metisn",options.ordering)) {
	  ierr=ZHPDperm_metis_n(A, rowscale, colscale, p, invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: metisn performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("metise",options.ordering)) {
          ierr=ZHPDperm_metis_e(A, rowscale, colscale, p, invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: metise performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       else if (!strcmp("amd",options.ordering)) {
          ierr=ZHPDperm_amd(A, rowscale, colscale, p, invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: amd performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("rcm",options.ordering)) {
	  ierr=ZHPDperm_rcm(A, rowscale, colscale, p, invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  ierr=ZHPDperm_null(A, rowscale, colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: no maximum weight matching, just scaling\n");fflush(stdout);
#endif
       }
       /* the scalings are real-valued anyway */
       for (i=0; i<n; i++) {
           pcolscale[i]=colscale[i].r;
       } /* end for */
       free(colscale);
 
       /* switch back to C to style */
       for (i=0; i<=n; i++)
	   A.ia[i]--;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]--;
       for (i=0; i<n; i++) {
	   p[i]--;
	   invq[i]--;
       } /* end for i */
       
#ifdef PRINT_CHECK
       /*
       for (j=0; j<n; j++)
           mexPrintf("%4ld",p[j]);
       mexPrintf("\n");
       for (j=0; j<n; j++)
           mexPrintf("%4ld",invq[j]);
       mexPrintf("\n");
       */
       for (j=0; j<n; j++) {
           if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n)
	      mexPrintf("ZHPDbilu ERROR, wrong p/invq after initial scaling, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j], invq[p[j]]);
       }  
#endif

#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
    
#ifdef PRINT_INFO2
       mexPrintf("p\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8ld",p[i]);
       }
       mexPrintf("\n");
       mexPrintf("invq (possibly extended, C-style)\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8ld",invq[i]);
       }
       mexPrintf("\n");
       mexPrintf("scaling that will be used\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8.1le",pcolscale[i]);
       }
       mexPrintf("\n");
#endif

    
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: matrix for BILU prepared\n");fflush(stdout);
#endif
#ifdef PRINT_INFO2
       mexPrintf("p\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8d",p[i]);
       }
       mexPrintf("\n");
       mexPrintf("invq\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8d",invq[i]);
       }
       mexPrintf("\n");
       mexPrintf("scaling\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8.1e",pcolscale[i]);
       }
       mexPrintf("\n");
#endif
#ifdef PRINT_INFO2
       mexPrintf("matrix C\n");
       for (i=0; i<C.nr; i++) {
	   mexPrintf("row %d\n", i+1);
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8d",C.rowind[i][j]);
	   }
	   mexPrintf("\n");
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8.1le",C.val[i][j].r);
           }
	   mexPrintf("\n");
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8.1le",C.val[i][j].i);
	   }
	   mexPrintf("\n");
       }
#endif
    }
    else { /* cosine successful or user-defined blocking */
       /* switch from C to FORTRAN style to apply ILUPACK scaling driver */
       for (i=0; i<=n; i++)
	   A.ia[i]++;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]++;
       /* permutation, right now ignored, but later used for companion matrix */
       p2   =(integer *)malloc((size_t)n*sizeof(integer));
       invq2=(integer *)malloc((size_t)n*sizeof(integer));
       ierr=ZHPDperm_null(A, rowscale,colscale, p2,invq2, &n, &options);
       /* the scalings are real-valued anyway */
       for (i=0; i<n; i++) {
           pcolscale[i]=colscale[i].r;
       } /* end for */
       /* switch back to C to style */
       for (i=0; i<=n; i++)
	   A.ia[i]--;
       for (i=0; i<A.ia[n]; i++)
	   A.ja[i]--;
       
       /* compute a companion matrix B of A(p,p) compressing the associated blocks */
       invblock=(integer *)malloc((size_t)n*sizeof(integer));
       j=0;
       for (i=0; i<mynblocks; i++) {
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         invblock[j]=i;
		 j++;
	   } /* end while */
       } /* end for i */

       B.nr=B.nc=mynblocks;
       B.ia=(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
       idx   =(integer *)malloc((size_t)mynblocks*sizeof(integer));
       idxpos=(integer *)calloc((size_t)mynblocks,sizeof(integer));
       cnt=0;
       
       /* first pass: detect nnz column by column */
       j=0;
       nnz=0;
       for (i=0; i<mynblocks; i++) {
	   /* scan block column i */
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         jj=p[j]; /* shortcut */
		 for (k=A.ia[jj]; k<A.ia[jj+1]; k++) {
		     ll=A.ja[k];
		     /* we will have q[l]=ll, so in total we are considering A(q[l],p[j]) */
		     l=invq[ll];
		     m=invblock[l];
		     if (!idxpos[m]) {
		        idx[cnt]=m;
			idxpos[m]=++cnt;
		     } /* end if */
		 } /* end for k */
		 j++;
	   } /* end while */
	   B.ia[i+1]=cnt;
	   nnz+=cnt;
	   for (k=0; k<cnt; k++) {
	       m=idx[k];
	       idxpos[m]=0;
	   } /* end for k */
	   cnt=0;  
       } /* end for i */

       /* change B.ia from nnz to pointer */
       B.ia[0]=0;
       for (i=0; i<mynblocks; i++) 
	   B.ia[i+1]+=B.ia[i];
       /* FORTRAN notation for companion matrix B for ILUPACK reordering routines */
       for (i=0; i<=mynblocks; i++) 
	   B.ia[i]++;
       
       B.ja=(integer *)      malloc((size_t)nnz*sizeof(integer));
       B.a =(doublecomplex *)malloc((size_t)nnz*sizeof(doublecomplex));
       B.nnz=nnz;
       /* second pass: detect pattern column by column */
       j=0;
       for (i=0; i<mynblocks; i++) {
	   /* scan block column i */
	   bi=myblocksize[i];
	   r=j; s=r+bi;
	   while (j<s) {
	         jj=p[j]; /* shortcut */
		 for (k=A.ia[jj]; k<A.ia[jj+1]; k++) {
		     ll=A.ja[k];
		     /* we will have q[l]=ll, so in total we are considering A(q[l],p[j]) */
		     l=invq[ll];
		     m=invblock[l];
		     if (!idxpos[m]) {
		        idx[cnt]=m;
			idxpos[m]=++cnt;
		     } /* end if */
		 } /* end for k */
		 j++;
	   } /* end while */
	   /* remember FORTRAN notation for B */
	   l=B.ia[i]-1;
	   for (k=0; k<cnt; k++) {
	       m=idx[k];
	       idxpos[m]=0;
	       B.ja[l] =m+1;
	       B.a[l].r  =1.0; /* dummy */
	       B.a[l++].i=0.0; /* dummy */
	   } /* end for k */
	   cnt=0;  
       } /* end for i */  
       free(idx);
       free(idxpos);
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: companion matrix computed\n");fflush(stdout);
#endif

#ifdef PRINT_INFO2
       mexPrintf("ZHPDbilu: call (scaling+) reordering\n");fflush(stdout);
       for (i=0; i<B.nc; i++){
	   mexPrintf("column %8ld\n",i+1);
	   for (j=B.ia[i]; j<B.ia[i+1]; j++){
	       k=B.ja[j-1];
	       mexPrintf("%8ld",k);
	   }
	   mexPrintf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++){
	       mexPrintf("%8.1le",B.a[j-1].r);
	   }
	   mexPrintf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++){
	       mexPrintf("%8.1le",B.a[j-1].i);
	   }
	   mexPrintf("\n");
       }
#endif

       /* select permutation and scaling driver by user request */
       nB=mynblocks;
       if (!strcmp("mtmetis",options.ordering)) {
          ierr=ZHPDperm_mtmetis(B, rowscale, colscale, p2, invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: mtmetis performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       else if (!strcmp("metisn",options.ordering)) {
	  ierr=ZHPDperm_metis_n(B, rowscale, colscale, p2, invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: metisn performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("metise",options.ordering)) {
          ierr=ZHPDperm_metis_e(B, rowscale, colscale, p2, invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: metise performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       else if (!strcmp("amd",options.ordering)) {
          ierr=ZHPDperm_amd(B, rowscale, colscale, p2, invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: amd performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("rcm",options.ordering)) {
	  ierr=ZHPDperm_rcm(B, rowscale, colscale, p2, invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  for (i=0; i<nB; i++)
	      p2[i]=invq2[i]=i+1;
#ifdef PRINT_INFO
	  mexPrintf("ZHPDbilu: no symmetric reordering used\n");fflush(stdout);
#endif
       }
       /* complex-valued scalings are not needed anymore */
       free(colscale);
       
       /* B not required anymore */
       free(B.ia);
       free(B.ja);
       free(B.a);
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
    
#ifdef PRINT_INFO2
       mexPrintf("block p as returned\n");
       for (i=0; i<B.nc; i++) {
	   mexPrintf("%8ld",p2[i]);
       }
       mexPrintf("\n");
       fflush(stdout);
       mexPrintf("block invq as returned \n");
       for (i=0; i<B.nc; i++) {
	   mexPrintf("%8ld",invq2[i]);
       }
       mexPrintf("\n");
       fflush(stdout);
#endif



       /* expand permutations from block to scalar */
       invq3=(integer *)malloc((size_t)n*sizeof(integer));
       p3   =(integer *)malloc((size_t)n*sizeof(integer));
       blockind=(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
       /* extract scalar starting positions from blocksize */
       blockind[0]=0;
       for (i=0; i<mynblocks; i++)
	   blockind[i+1]=blockind[i]+myblocksize[i];
       /* expand block p2 to scalar p3 */
       j=0;
       for (i=0; i<mynblocks; i++) {
	   k=p2[i]-1; /* FORTRAN to C notation */
	   for (l=blockind[k]; l<blockind[k+1]; l++)
	       p3[j++]=l;
       } /* end for i */

    
       /* blocksize is now incoroprated by the permutations */
       if (options.blocksize==NULL && !ilu1t) { 
	  /* no further a priori block requested, transfer block information */
	  /* use invblock as buffer */
	  for (i=0; i<mynblocks; i++) {
	      k=p2[i]-1; /* FORTRAN to C notation */
	      invblock[i]=myblocksize[k];
	  }
	  memcpy(myblocksize, invblock, mynblocks*sizeof(integer));
       }
#ifdef PRINT_CHECK
       if (mynblocks>0 && myblocksize!=NULL) {
	  k=0;
	  for (j=0; j<mynblocks; j++) {
	      k+=myblocksize[j];
	      if (myblocksize[j]<=0 || myblocksize[j]>n)
		 printf("ZHPDbilu ERROR, wrong block size, permuted cosine, step %ld, blocksize[%ld]=%ld\n", j,j,myblocksize[j]);
	      if (k>n)
		 printf("ZHPDbilu ERROR, cummulative block size exceeds system size, permuted cosine, step %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
	  }  
	  if (k!=n)
	     printf("ZHPDbilu ERROR, cummulative block does not match system size, permuted cosine, nblocks %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
	  for (j=0; j<mynblocks; j++)
	      printf("%4ld",myblocksize[j]);
	  printf("\n");fflush(stdout);
       }
#endif

    
       /* compute block q2 from block invq2 and use p2 as buffer */
       for (i=0; i<mynblocks; i++)
	   p2[invq2[i]-1]=i;
       j=0;
       for (i=0; i<mynblocks; i++) {
	   k=p2[i]; /* now q2 (stored in p2) already uses C notation */
	   for (l=blockind[k]; l<blockind[k+1]; l++)
	       invblock[j++]=l; /* compute scalar q3 and use invblock as buffer */
       } /* end for i */
       /* compute scalar invq3 from scalar q3 (stored in invblock) */
       for (i=0; i<n; i++)
	   invq3[invblock[i]]=i;

#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: permutation expanded\n");fflush(stdout);
#endif
    
    
       /* block structures are not needed anymore */
       free(blockind);
       free(invq2);
       free(p2);
    
       /* replace permutations by their product with the symmetric permutations */
       /* p <- p[p3] use invblock as buffer */
       for (i=0; i<n; i++)
	   invblock[i]=p[p3[i]];
       memcpy(p, invblock, n*sizeof(integer));
       /* invq <- invq3[invq] use invblock as buffer */
       for (i=0; i<n; i++)
           invblock[i]=invq3[invq[i]];
       memcpy(invq, invblock, n*sizeof(integer));
#ifdef PRINT_CHECK
       for (j=0; j<n; j++)
	   mexPrintf("%4ld",p[j]);
       mexPrintf("\n");
       for (j=0; j<n; j++)
	   mexPrintf("%4ld",invq[j]);
       mexPrintf("\n");
       for (j=0; j<n; j++) {
	   if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n)
	      mexPrintf("ZHPDbilu ERROR, wrong p/invq after reorering the compressed graph, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j],invq[p[j]]);
       }  
#endif
    
       free(invblock);
       free(p3);
       free(invq3);
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: permutation product built\n");fflush(stdout);
#endif
#ifdef PRINT_INFO2
       mexPrintf("p (possibly extended, C-style)\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8ld",p[i]);
       }
       mexPrintf("\n");
       mexPrintf("invq (possibly extended, C-style)\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8ld",invq[i]);
       }
       mexPrintf("\n");
       mexPrintf("scaling that will be used\n");
       for (i=0; i<n; i++) {
	   mexPrintf("%8.1le",pcolscale[i]);
       }
       mexPrintf("\n");
#endif

    
#ifdef PRINT_INFO
       mexPrintf("ZHPDbilu: matrix for BILU prepared\n");fflush(stdout);
#endif
#ifdef PRINT_INFO2
       mexPrintf("p\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8d",p[i]);
       }
       mexPrintf("\n");
       mexPrintf("invq\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8d",invq[i]);
       }
       mexPrintf("\n");
       mexPrintf("scaling\n");
       for (i=0; i<A.nr; i++) {
	   mexPrintf("%8.1e",pcolscale[i]);
       }
       mexPrintf("\n");
#endif
#ifdef PRINT_INFO2
       mexPrintf("matrix C\n");
       for (i=0; i<C.nr; i++) {
	   mexPrintf("row %d\n", i+1);
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8d",C.rowind[i][j]);
	   }
	   mexPrintf("\n");
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8.1le",C.val[i][j].r);
	   }
	   mexPrintf("\n");
	   for (j=0; j<C.ncol[i]; j++) {
	       mexPrintf("%8.1le",C.val[i][j].i);
	   }
	   mexPrintf("\n");
       }
#endif
    } /* end else */
    

    /* start already exporting the scalings */
    /* 6. lhs: S_L */
    plhs[3]=mxCreateSparse((mwSize)n,(mwSize)n, (mwSize)n, mxREAL);
    S_L=plhs[3];
    ja     =(mwIndex *)mxGetIr(S_L);
    ia     =(mwIndex *)mxGetJc(S_L);
    valuesR=(double *) mxGetPr(S_L);
    for (i=0; i<n; i++) {
        ia[i]=i;
        ja[i]=i;
	valuesR[i]=pcolscale[i];
    } /* end for i */
    ia[n]=n;
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: S_L already exported\n");fflush(stdout);
#endif

    /* pcolscale not required anymore since C is already properly scaled */
    free(pcolscale);
    
#ifdef PRINT_INFO
    mexPrintf("call BILUC, droptol=%8.1le\n", options.droptol); fflush(stdout);
#endif


#ifdef PRINT_INFO2
    mexPrintf("matrix C to be passed\n");
    for (i=0; i<C.nc; i++) {
        mexPrintf("column %ld\n", i);
        for (j=0; j<C.ncol[i]; j++) {
	    mexPrintf("%16ld",C.rowind[i][j]);
	}
        mexPrintf("\n");
        for (j=0; j<C.ncol[i]; j++) {
	    mexPrintf("%16.8le",C.val[i][j].r);
	}
        mexPrintf("\n");
        for (j=0; j<C.ncol[i]; j++) {
	    mexPrintf("%16.8le",C.val[i][j].i);
	}
        mexPrintf("\n");
    }
    mexPrintf("p\n"); 
    for (i=0; i<C.nc; i++)
        mexPrintf("%16ld",p[i]);
    mexPrintf("\n");
    mexPrintf("invq\n"); 
    for (i=0; i<C.nc; i++)
        mexPrintf("%16ld",invq[i]);
    mexPrintf("\n");
#endif

    if (options.blocksize==NULL && !cosine) {
       if (myblocksize!=NULL)
	  free(myblocksize);
       myblocksize=NULL;
       mynblocks=0;
    }					   
    if (options.blocksize==NULL) { 
       if (ilu1t) {
	  if (myblocksize!=NULL)
	     free(myblocksize);
	  myblocksize=NULL;
	  mynblocks=0;
       }
    }
    else {
       myblocksize=options.blocksize;
       mynblocks=options.nblocks;
    }

    ierr=ZHPDbilu(&C, &BL,&BiD,&BL, NULL,NULL, p,invq, &determinant, myblocksize,mynblocks,
		  options.droptol,ilu1t,pa,perturbation,invert_blocks,level_of_fill);
    for (i=0; i<n; i++)
        determinant.r-=2*log(valuesR[i]);

#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: bilu completed, ierr=%d\n",ierr);fflush(stdout);
#endif
    if (options.blocksize==NULL && !ilu1t) {
       if (myblocksize!=NULL)
	  free(myblocksize);
       myblocksize=NULL;
       mynblocks=0;
    }
    /* matrix C and A are not longer required */
    free(A.ia);  
    free(A.ja);
    free(A.a);
    free(C.ncol);  
    free(C.rowind);
    free(C.val);

    if (ierr) {
       mexErrMsgTxt("ZHPDbilu: block ILU abnormally terminated") ;
    }
    /* 4. lhs: P */
    plhs[2]=mxCreateSparse((mwSize)n,(mwSize)n, (mwSize)n, mxREAL);
    P=plhs[2];
    ja     =(mwIndex *)mxGetIr(P);
    ia     =(mwIndex *)mxGetJc(P);
    valuesR=(double *) mxGetPr(P);
    for (i=0; i<n; i++) {
        ia[i]=i;
        ja[i]=p[i];
	valuesR[i]=1.0;
    } /* end for i */
    ia[n]=n;
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: P exported\n");fflush(stdout);
#endif
    free(invq);
    free(p);

    
    /* 4. lhs: determinant */
    plhs[4]=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxCOMPLEX);
    pr=mxGetPr(plhs[4]);
    *pr=determinant.r;
    pr=mxGetPi(plhs[4]);
    *pr=determinant.i;
    
    
    /* export options */
    if (options.blocksize!=NULL) {
       plhs[5]=mxCreateStructMatrix((mwSize)1, (mwSize)1, 9, optionsnamesx);
    }
    else {
       plhs[5]=mxCreateStructMatrix((mwSize)1, (mwSize)1, 8, optionsnames);
    }
    options_output=plhs[5];
    /* field 0: cosine */
    if (cosine) {
       output_buf=(char *)mxCalloc((size_t)4, (size_t)sizeof(char));
       strcpy(output_buf,"yes");
    }
    else {
       output_buf=(char *)mxCalloc((size_t)3, (size_t)sizeof(char));
       strcpy(output_buf,"no");
    }
    fout = mxCreateString(output_buf);
    mxSetFieldByNumber(options_output, (mwSize)0, 0, fout);
    mxFree(output_buf);
    /* field 1: ilu1t */
    if (ilu1t==BLOCK_ILU1T) {
       output_buf=(char *)mxCalloc((size_t)6, (size_t)sizeof(char));
       strcpy(output_buf,"ilu1t");
    }
    else if (ilu1t==BLOCK_ILUPT) {
       output_buf=(char *)mxCalloc((size_t)6, (size_t)sizeof(char));
       strcpy(output_buf,"ilupt");
    }
    else if (ilu1t==BLOCK_SUPERNODES) {
       output_buf=(char *)mxCalloc((size_t)11, (size_t)sizeof(char));
       strcpy(output_buf,"supernodes");
    }
    else if (ilu1t==BLOCK_APP_SUPERNODES) {
       output_buf=(char *)mxCalloc((size_t)9, (size_t)sizeof(char));
       strcpy(output_buf,"appsuper");
    }
    else {
       output_buf=(char *)mxCalloc((size_t)5, (size_t)sizeof(char));
       strcpy(output_buf,"none");
    }
    fout = mxCreateString(output_buf);
    mxSetFieldByNumber(options_output, (mwSize)0, 1, fout);
    mxFree(output_buf);
    /* field 2: progressive_aggregation */
    if (pa) {
       output_buf=(char *)mxCalloc((size_t)4, (size_t)sizeof(char));
       strcpy(output_buf,"yes");
    }
    else {
       output_buf=(char *)mxCalloc((size_t)3, (size_t)sizeof(char));
       strcpy(output_buf,"no");
    }
    fout = mxCreateString(output_buf);
    mxSetFieldByNumber(options_output, (mwSize)0, 2, fout);
    mxFree(output_buf);
    /* field 4: droptol */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=options.droptol;
    mxSetFieldByNumber(options_output, (mwSize)0, 4, tmp);
    /* field 3: ordering */    
    output_buf=(char *)mxCalloc((size_t)strlen(options.ordering)+1, (size_t)sizeof(char));
    strcpy(output_buf,options.ordering);
    if (flag)
       free(options.ordering);
    fout = mxCreateString(output_buf);
    mxSetFieldByNumber(options_output, (mwSize)0, 3, fout);
    mxFree(output_buf);
    /* field 5: perturbation */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=perturbation;
    mxSetFieldByNumber(options_output, (mwSize)0, 5, tmp);
    /* field 6: invert_blocks */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=invert_blocks;
    mxSetFieldByNumber(options_output, (mwSize)0, 6, tmp);
    /* field 7: level_of_fill */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=level_of_fill;
    mxSetFieldByNumber(options_output, (mwSize)0, 7, tmp);
    /* field 8: blocksize (optionally) */    
    if (options.blocksize!=NULL) {
       tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)options.nblocks, mxREAL);
       pr=mxGetPr(tmp);
       for (i=0; i<options.nblocks; i++) {
	   pr[i]=options.blocksize[i];
       } /* end for i */
       free(options.blocksize);
       mxSetFieldByNumber(options_output, (mwSize)0, 8, tmp);
    }
    
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: options structure exported\n");fflush(stdout);
#endif


    /* 1.lhs: export BL */
    mydims[0]=BL.nblocks;
    plhs[0]=mxCreateCellArray((mwSize)1, mydims);
    BL_output=plhs[0];
    for (i=0; i<BL.nblocks; i++) {
        /* set up new block column for BL{i} with four elements J, I, L, D */
        Block=mxCreateStructMatrix((mwSize)1, (mwSize)1, 4, BLnames);

	/* structure element 0:  J */
	/* create BL{i}.J */
	m=BL.nblockcol[i];
	BlockJ=mxCreateDoubleMatrix((mwSize)1,(mwSize)m, mxREAL);
	/* copy data */
	pr=(double *)mxGetPr(BlockJ);
	for (j=0; j<m; j++)
	    pr[j]=BL.colind[i][j]+1;
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 0, BlockJ);

	/* structure element 1:  I */
	/* create empty BL{i}.I */
	k=BL.nblockrow[i];
	BlockI=mxCreateDoubleMatrix((mwSize)1,(mwSize)k, mxREAL);
	/* copy data */
	pr=(double *)mxGetPr(BlockI);
	for (j=0; j<k; j++)
	    pr[j]=BL.rowind[i][j]+1;
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 1, BlockI);

	/* structure element 2:  L */
	/* create empty BL{i}.L */
	BlockL=mxCreateDoubleMatrix((mwSize)k,(mwSize)m, mxCOMPLEX);
	pr=(double *)mxGetPr(BlockL);
	pi=(double *)mxGetPi(BlockL);
	for (j=0; j<k*m; j++) {
	    pr[j]=BL.valE[i][j].r;
	    pi[j]=BL.valE[i][j].i;
	}
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 2, BlockL);

	/* structure element 3:  D */
	/* create identity sparse m x m matrix BL{i}.D */
	BlockD=mxCreateSparse((mwSize)m,(mwSize)m, (mwSize)m, mxREAL);
	ja     =(mwIndex *)mxGetIr(BlockD);
	ia     =(mwIndex *)mxGetJc(BlockD);
	valuesR=(double *) mxGetPr(BlockD);
	for (j=0; j<m; j++) {
	    ia[j]=j;
	    ja[j]=j;
	    valuesR[j]=1.0;
	} /* end for j */
	ia[m]=m;
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 3, BlockD);

	/* finally set output BL{i} */
	mxSetCell(BL_output, (mwIndex)i,Block);
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: BL exported\n");fflush(stdout);
#endif
    

    

    /* export BiD */
    mydims[0]=BiD.nblocks;
    plhs[1]=mxCreateCellArray((mwSize)1, mydims);
    BiD_output=plhs[1];
    for (i=0; i<BiD.nblocks; i++) {
        /* set up new block column for BiD{i} with four elements J, I, L, D */
        Block=mxCreateStructMatrix((mwSize)1, (mwSize)1, 4, BLnames);

	/* structure element 0:  J */
	/* create BiD{i}.J */
	m=BiD.nblockcol[i];
	BlockJ=mxCreateDoubleMatrix((mwSize)1,(mwSize)m, mxREAL);
	/* copy data */
	pr=(double *)mxGetPr(BlockJ);
	for (j=0; j<m; j++)
	    pr[j]=BiD.colind[i][j]+1;
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 0, BlockJ);

	/* structure element 1:  I */
	/* create empty BiD{i}.I */
	BlockI=mxCreateDoubleMatrix((mwSize)1,(mwSize)0, mxREAL);
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 1, BlockI);

	/* structure element 2:  L */
	/* create empty BiD{i}.L */
	BlockL=mxCreateDoubleMatrix((mwSize)0,(mwSize)m, mxREAL);
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 2, BlockL);

	/* structure element 3:  D */
	/* create dense m x m matrix BiD{i}.D */
	BlockD=mxCreateDoubleMatrix((mwSize)m,(mwSize)m, mxCOMPLEX);
	pr=(double *)mxGetPr(BlockD);
	pi=(double *)mxGetPi(BlockD);
	for (j=0; j<m*m; j++) {
	    pr[j]=BiD.valD[i][j].r;
	    pi[j]=BiD.valD[i][j].i;
	}
	if (!invert_blocks) {
	   /* scalar case: pass psychological "Cholesky factor" sqrt(dii) rather than 1/dii */
	   if (m==1)
	      pr[0]=1/sqrt(pr[0]);
	   else { /* set strict upper triangular part to 0 */
	      pr+=m;
	      pi+=m;
	      for (j=1; j<m; j++) {
		  for (k=0; k<j; k++)
		      pr[k]=pi[k]=0.0;
		  pr+=m;
		  pi+=m;
	      } /* end for k */
	   } /* end if-else m=1 */
	} /* end if */
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 3, BlockD);

	/* finally set output BiD{i} */
	mxSetCell(BiD_output, (mwIndex)i,Block);
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: BiD exported\n");fflush(stdout);
#endif

    /* finally release auxiliary memory */
    ZSparseBlockDelete(&BL);
    ZSparseBlockDelete(&BiD);


#ifdef PRINT_INFO
    mexPrintf("ZHPDbilu: memory released\n");fflush(stdout);
#endif
    return;       

}
