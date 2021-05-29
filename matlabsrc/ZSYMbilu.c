/* $Id: ZSYMbilu.c 7315 2021-05-28 21:00:20Z bolle $ */
/* ========================================================================== */
/* === ZSYMbilu mexFunction ================================================= */
/* ========================================================================== */

/*
    Usage:

    Given a permutation matrix P as well as a scaling matrix S_L,
    compute a block-structured approximate factorization

          P^T S_L A S_L P ~  BL BiD^{-1} BL^T

    with a block (inverse) diagonal matrix BiD, block unit lower triangular 
    factor BL.


    Example:

    [BL,BiD,P,S_L,determinant,pivots,options]=ZSYMbilu(A,options)


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	March 14, 2017. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.

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
#define USE_COSINE_BLOCKS /*  */
#define TAU      0.8
/* #define PRINT_CHECK */
/* #define PRINT_INFO  */
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif

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
                       *optionsnames[] ={"cosine","blocking_strategy","progressive_aggregation","matching","ordering","droptol","perturbation","isdefinite","invert_blocks","level_of_fill"},
		       *optionsnamesx[]={"cosine","blocking_strategy","progressive_aggregation","matching","ordering","droptol","perturbation","isdefinite","invert_blocks","level_of_fill","blocksize"};
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
    integer            i,j,k,l,ii,ll,jj,m,n,*p,*invq,*invq2,*p2,*invq3,*p3,r,s,mynblocks,*myblocksize,flag=0,*pivots,
                       nB=0, ierr,bi, *invblock=NULL, *blockind, cosine=0, ilu1t=1, pa=1,
                       perturbation=1, isdefinite=0,
                       invert_blocks=1,
                       level_of_fill=1,
                       cnt, *idx, *idxpos;
    double             *valuesR, *pr, *valuesI, *pi,
                       *pcolscale, *prowscale, val,
                       *pcolscale2, *prowscale2;
    doublecomplex      *colscale, *rowscale;
    doublecomplex      determinant;
    size_t             mrows, ncols;
    ZILUPACKparam      options;
    Zmat               A,B;
    ZSparseMatrix      C;
    ZSparseBlockMatrix BL,BiD;
#ifdef PRINT_CHECK
    integer *check_buff;
#endif
#ifdef _PROFILING_
  double timeBegin,
         time_total=0,
         time_bilu=0,
         time_pre_processing=0,
         time_cosine=0,
         time_reordering=0,
         time_post_processing=0;

  timeBegin=time_total=omp_get_wtime();
#endif
    
    if (nrhs!=2)
       mexErrMsgTxt("Two input arguments required.");
    else if (nlhs!=7)
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
       mexErrMsgTxt("ZSYMbilu: input matrix must be in sparse format.") ;
    }
    n=mrows;
#ifdef PRINT_CHECK
    check_buff=(integer *)calloc(n,sizeof(integer));;
#endif

    /* copy input matrix to sparse row format */
    A.nc=A.nr=n;
    A.ia=(integer *)      malloc((size_t)(n+1)*sizeof(integer));
    A.ja=(integer *)      malloc((size_t)nnz     *sizeof(integer));
    A. a=(doublecomplex *)malloc((size_t)nnz     *sizeof(doublecomplex));

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

#ifdef PRINT_CHECK0
    for (i=0; i<A.nr; i++) {
        mexPrintf("row %6ld:\n",i);
	for (k=A.ia[i]; k<A.ia[i+1]; k++)
	    mexPrintf("%8ld",A.ja[k-1]-1);
        mexPrintf("\n");
	for (k=A.ia[i]; k<A.ia[i+1]; k++)
	    mexPrintf("%8.1le",A.a[k-1].r);
	for (k=A.ia[i]; k<A.ia[i+1]; k++)
	    mexPrintf("%8.1le",A.a[k-1].i);
        mexPrintf("\n");
    }
#endif
    
    /* initialize ILUPACK options structure to its default options */
    ZSYMAMGinit(&A,&options);
#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: sparse matrix imported\n");fflush(stdout);
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

    options.nblocks=0;
    options.blocksize=NULL;
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
	   if (!strcmp("matching",fnames[ifield])) {
	      options.matching=*mxGetPr(tmp);
	   }
	   else if (!strcmp("droptol",fnames[ifield])) {
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
		 mexErrMsgTxt("ZSYMbilu: options.blocksize must exactly match the size of the matrix") ;
	      }
	      cnt=0;
	   }
	   else {
	      /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	}
    }
#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: options structure imported\n");fflush(stdout);
    mexPrintf("ZSYMbilu: cosine=%ld, ilu1t=%ld, progressive_aggregation=%ld\n",
	      cosine,ilu1t,pa);fflush(stdout);
#endif
    mxFree(classIDflags);
    mxFree(fnames);

    
    /* memory for  permutation */
    p   =(integer *)malloc((size_t)n*sizeof(integer));
    invq=(integer *)malloc((size_t)n*sizeof(integer));
    /* row scaling vector */
    pcolscale=(double *)malloc((size_t)n*sizeof(double));
    prowscale=pcolscale;
    colscale=(doublecomplex *)malloc((size_t)n*sizeof(doublecomplex));
    rowscale=colscale;
    nB=n;
#ifdef _PROFILING_
    time_pre_processing+=omp_get_wtime()-timeBegin;
#endif



    
    
    if (options.nblocks==n||options.nblocks==0||options.blocksize==NULL) {
       /* if matching is requested */
       if (options.matching && cosine) {
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* call maximum weight matching followed by a symmetrizing routine */
	  ierr=ZSYMperm_mc64_null(A, rowscale,colscale, p,invq, &nB, &options);
	  /* now the permuted system A(p,p) refers to a matrix, where the leading
	     nB x nB block refers to 1x1 blocks whereas the remaining (n-nB) x (n-nB) 
	     prefers to work with 2x2 pivots. This block is to be maintained in the
	     sequel.
	     Note also that A has been rescaled by some D_rowscale*A*D_colscale
	  */
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
              pcolscale[i]=colscale[i].r;
	  } /* end for */
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  } /* end for i */
#ifdef PRINT_CHECK
	  for (i=0; i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;	      
	  for (i=0; i<n; i++) {
	      if (invq[p[i]]!=i) {
	         mexPrintf("ZSYMbilu ERROR, inverse permutation does not match, step %ld, p[%ld]=%ld, invq[p[%ld]]=%ld\n",i,i,p[i],i,invq[p[i]]);
		 fflush(stdout);
	      }
	  }
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
       
	  /* if we want to use the cosine blocking strategy on top of matching,
	     we first have to replace A by a compressed counter part
	  */

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* switch matrix A back to C notation and adapt data structures */
	  for (i=0; i<=n; i++)
	      A.ia[i]--;
	  for (i=0; i<A.ia[n]; i++)
	      A.ja[i]--;
	  /* start by copying the permuted pattern */
	  C.nc=C.nr=n;
	  C.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
	  C.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
	  C.val   =(doublecomplex **) malloc((size_t)n*sizeof(doublecomplex *));
	     
	  for (j=0; j<n; j++) {
	      ii=p[j];
	      C.ncol[j]=A.ia[ii+1]-A.ia[ii];
	      C.rowind[j]=(integer *)malloc((size_t)C.ncol[j]*sizeof(integer));
	      C.val[j]   =(doublecomplex *) malloc((size_t)C.ncol[j]*sizeof(doublecomplex));
	      for (k=0,i=A.ia[ii]; k<C.ncol[j]; k++,i++) {
		  C.rowind[j][k]=invq[A.ja[i]];
	      }
	  } /* end for j */
#ifdef PRINT_INFO
	  /* DPrintMatrix(C); */
	  mexPrintf("ZSYMbilu: permuted sparse matrix data structures set up\n");fflush(stdout);
#endif


	  /* compress matrix, use invq temporarily as buffer */
	  /* the leading nB rows correspond to 1x1 pivots */
	  for (i=0; i<nB; i++)
	      invq[i]=i;
	  j=nB;
	  for (; i<n; i+=2)
	      invq[i]=invq[i+1]=j++;

	  /* overwrite row indices by their compressed representation, possibly leading to
	     duplicate indices 
	  */
	  for (j=0; j<n; j++) {
	      for (k=0; k<C.ncol[j]; k++) {
		  i=C.rowind[j][k];
		  C.rowind[j][k]=invq[i];
	      } /* end for k */
	  } /* end for j */

	  /* use invq temporarily as buffer for flags */
	  for (i=0; i<n; i++)
	      invq[i]=0;

	  /* skip duplicate entries in the leading nB rows */
	  for (j=0; j<nB; j++) {
	      /* new start of column j */
	      k=0;
	      for (i=0; i<C.ncol[j]; i++) {
		  ii=C.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     C.rowind[j][k++]=ii;
		  } /* end if */
	      } /* end for i */
	      /* clear check mark array */
	      for (i=0; i<k; i++) {
		  ii=C.rowind[j][i];
		  invq[ii]=0;
	      }
	      /* update number of nonzero row indices in column j */
	      C.ncol[j]=k;
	      /* reduce amount of memory */
	      C.rowind[j]=(integer *)      realloc(C.rowind[j],(size_t)k*sizeof(integer));
	      C.val[j]   =(doublecomplex *)realloc(C.val[j],   (size_t)k*sizeof(doublecomplex));
	  } /* end for j */
	  
	  /* merge duplicate entries and rows in the remaining n-nB rows */
	  for (m=nB; j<n; m++,j++) {
	      /* new start of column m */
	      k=0;

	      /* make sure that the current column is long enough */
	      C.rowind[m]=(integer *)realloc(C.rowind[m],(size_t)(C.ncol[j]+C.ncol[j+1])*sizeof(integer));
	      /* first scan column j */
	      for (i=0; i<C.ncol[j]; i++) {
		  ii=C.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     C.rowind[m][k++]=ii;
		  } /* end if */
	      } /* end for i */

	      j++;
	      /* second scan column j+1 */
	      for (i=0; i<C.ncol[j]; i++) {
		  ii=C.rowind[j][i];
		  /* entry is not duplicate yet */
		  if (!invq[ii]) {
		     /* check mark */
		     invq[ii]=-1;
		     /* shift entry */
		     C.rowind[m][k++]=ii;
		  } /* end if */
	      } /* end for i */
	     
	      /* clear check mark array of column j and j+1 */
	      for (i=0; i<k; i++) {
		  ii=C.rowind[m][i];
		  invq[ii]=0;
	      }
	      /* update number of nonzero row indices in column j */
	      C.ncol[m]=k;
	      /* reduce amount of memory */
	      C.rowind[m]=(integer *)      realloc(C.rowind[m],(size_t)k*sizeof(integer));
	      C.val[m]   =(doublecomplex *)realloc(C.val[m],   (size_t)k*sizeof(doublecomplex));
	  } /* end for m */
	  C.nr=C.nc=m;

	  /* get rid of the remaining columns */
	  for (j=m; j<n; j++) {
	      free(C.rowind[j]);
	      free(C.val[j]);
	  }
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: companion matrix for cosine-based algorithm computed\n");fflush(stdout);
#endif

	  /* apply the cosine-based blocking to the compressed matrix */
	  myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
	  p2         =(integer *)malloc((size_t)n*sizeof(integer));
	  mynblocks=0;
	  for (j=0; j<m; j++)
	      p2[j]=invq[j]=j;
	  Zmcosine_sblocks(&C, p2,invq, myblocksize, &mynblocks, TAU);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: cosine applied\n");fflush(stdout);
#endif
#ifdef PRINT_CHECK
	  for (i=0; i<m; i++) {
	      if (p2[i]<0 || p2[i]>=m) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=m) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;
	  for (i=0; i<m; i++) {
	      if (check_buff[p2[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
		 fflush(stdout);
	      }
	      check_buff[p2[i]]=-1;
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;
	  for (i=0; i<m; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<m; i++)
	      check_buff[i]=0;	      
	  for (i=0; i<m; i++) {
	      if (invq[p2[i]]!=i) {
	         mexPrintf("ZSYMbilu ERROR, inverse permutation does not match, step %ld, p2[%ld]=%ld, invq[p2[%ld]]=%ld\n",i,i,p2[i],i,invq[p2[i]]);
		 fflush(stdout);
	      }
	  }
#endif

	  /* prolongate permutation to its original size, use invq as auxiliary vector */
	  /* compute inverse mapping index->block */
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

	  /* expand blocks of B(p2,p2) */
	  j=0;
	  for (i=0;i<m; i++) {
	      /* consider B(:,k), B(k,:) */
	      k=p2[i];
	      /* 1x1 block, no expansion */
	      if (k<nB)
		 invq[j++]=k;
	      else { /* 2x2 block, k refers to a block of two rows/columns */
		 myblocksize[invblock[i]]++;
		 invq[j]=2*k-nB;
		 invq[j+1]=invq[j]+1;
		 j+=2;
	      }
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong expanded permutation, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
	  

	  /* build product permutation for the columns (symmetric maximum 
	     weighted matching followed cosine-based blocking)
	  */
	  for (i=0; i<n; i++) 
	      p2[i]=p[invq[i]];
	  /* rewrite permutation */
	  for (i=0; i<n; i++)
	      p[i]=p2[i];
	  /* associated inverse permutation */
	  for (i=0; i<n; i++)
	      invq[p[i]]=i;
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  /* remove remaining columns */
	  for (j=0; j<m; j++) {
	      free(C.rowind[j]);
	      free(C.val[j]);
	  }
	  free(p2);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: permutation after cosine expanded to incorporate 2x2 pivots\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_cosine+=omp_get_wtime()-timeBegin;
#endif

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* now set up the C for BILDL, recall that A is already shifted w.r.t. C-style,
	     further remember that we already allocated memory for C.ncol, C.rowind 
	  */
	  /* use slightly different data structure for bildl */
	  C.nr=A.nr; C.nc=n;
	  C.nnz=A.nnz;
	  for (i=0; i<n; i++) {
	      C.ncol[i]  =A.ia[i+1]-A.ia[i];
	      C.rowind[i]=A.ja+A.ia[i];
	      C.val   [i]=A.a +A.ia[i];
	  } /* end for i */
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: proper sparse matrix data structures based on A set up\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* now we have applied maximum weight matching followed by cosine-based blocking,
	     yet, a symmetric reordering taking into account the induced block structure
	     still has to be applied
	  */
       } /* end if matching and cosine */
       else if (options.matching) { /* only matching */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  if (!strcmp("mtmetis",options.ordering)) {
#ifdef _MC64_MATCHING_
	     ierr=ZSYMperm_mc64_mtmetis(A, rowscale,colscale, p,invq, &n, &options);
#else
	     ierr=ZSYMperm_matching_mtmetis(A, rowscale,colscale, p,invq, &n, &options);
#endif
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: matching+mtmetis performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  else if (!strcmp("metisn",options.ordering)) { /* METIS nested dissection by nodes */
#ifdef _MC64_MATCHING_
	     ierr=ZSYMperm_mc64_metis_n(A, rowscale,colscale, p,invq, &n, &options);
#else
	     ierr=ZSYMperm_matching_metis_n(A, rowscale,colscale, p,invq, &n, &options);
#endif	     
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: matching+metisn performed\n");fflush(stdout);
#endif
	  }
	  else if (!strcmp("metise",options.ordering)) {
#ifdef _MC64_MATCHING_
	     ierr=ZSYMperm_mc64_metis_e(A, rowscale,colscale, p,invq, &n, &options);
#else
	     ierr=ZSYMperm_matching_metis_e(A, rowscale,colscale, p,invq, &n, &options);
#endif
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: matching+metise performed\n");fflush(stdout);
#endif
	  }
#endif // _USE_METIS4_
	  else if (!strcmp("amd",options.ordering)) {
#ifdef _MC64_MATCHING_
	     ierr=ZSYMperm_mc64_amd(A, rowscale,colscale, p,invq, &n, &options);
#else
             ierr=ZSYMperm_matching_amd(A, rowscale,colscale, p,invq, &n, &options);
#endif
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: matching+amd performed\n");fflush(stdout);
#endif
	  }
	  else {
#ifdef _MC64_MATCHING_
	     ierr=ZSYMperm_mc64_rcm(A, rowscale,colscale, p,invq, &n, &options);
#else
	     ierr=ZSYMperm_matching_rcm(A, rowscale,colscale, p,invq, &n, &options);
#endif
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: matching+rcm performed\n");fflush(stdout);
#endif
	  }
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
	      pcolscale[i]=colscale[i].r;
	  } /* end for */
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  /* no cosine blocking */
	  myblocksize=NULL;
	  mynblocks=0;

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* set up the C for BILDL */
	  /* switch matrix A back to C notation and adapt data structures for bildl (and cosine_blocks) */
	  for (i=0; i<=n; i++)
	      A.ia[i]--;
	  for (i=0; i<A.ia[n]; i++)
	      A.ja[i]--;
	  /* use slightly different data structure for bildl */
	  C.nr=A.nr; C.nc=n;
	  C.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
	  C.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
	  C.val   =(doublecomplex **) malloc((size_t)n*sizeof(doublecomplex *));
	  C.nnz=A.nnz;
	  for (i=0; i<n; i++) {
	      C.ncol[i]  =A.ia[i+1]-A.ia[i];
	      C.rowind[i]=A.ja+A.ia[i];
	      C.val   [i]=A.a +A.ia[i];
	  } /* end for i */
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* we are ready to apply BILDL */
       } /* end if only matching */
       else if (cosine) { /* only cosine */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do some scaling in advance */
	  ierr=ZSYMperm_null(A, rowscale,colscale, p,invq, &n, &options);
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
	      pcolscale[i]=colscale[i].r;
	  } /* end for */
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif

	  /* apply the cosine-based blocking to the compressed matrix */
	  myblocksize=(integer *)malloc((size_t)n*sizeof(integer));
	  mynblocks=0;

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* set up the C for BILDL */
	  /* switch matrix A back to C notation and adapt data structures for bildl (and cosine_blocks) */
	  for (i=0; i<=n; i++)
	      A.ia[i]--;
	  for (i=0; i<A.ia[n]; i++)
	      A.ja[i]--;
	  /* use slightly different data structure for bildl */
	  C.nr=A.nr; C.nc=n;
	  C.ncol  =(integer *)       malloc((size_t)n*sizeof(integer));
	  C.rowind=(integer **)      malloc((size_t)n*sizeof(integer *));
	  C.val   =(doublecomplex **)malloc((size_t)n*sizeof(doublecomplex *));
	  C.nnz=A.nnz;
	  for (i=0; i<n; i++) {
	      C.ncol[i]  =A.ia[i+1]-A.ia[i];
	      C.rowind[i]=A.ja+A.ia[i];
	      C.val   [i]=A.a +A.ia[i];
	  } /* end for i */
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  Zmcosine_sblocks(&C, p,invq, myblocksize,&mynblocks, TAU);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: cosine applied\n");fflush(stdout);
#endif
#ifdef _PROFILING_
	  time_cosine+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
	  /* we have applied the cosine blocking only and still need to reorder
	     the system taking into account the associated block structure 
	  */
       } /* end if only cosine */
       else { /* neither cosine nor matching nor blocks */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do scaling and simple reordering without asking for matchings */
	  if (!strcmp("mtmetis",options.ordering)) {
	     ierr=ZSYMperm_mtmetis(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: mtmetis performed\n");fflush(stdout);
#endif
	  }
#ifdef _USE_METIS4_
	  else if (!strcmp("metisn",options.ordering)) {
	     ierr=ZSYMperm_metis_n(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: metisn performed\n");fflush(stdout);
#endif
	  }
	  else if (!strcmp("metise",options.ordering)) {
	     ierr=ZSYMperm_metis_e(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: metise performed\n");fflush(stdout);
#endif
	  }
#endif // _USE_METIS4_
	  else if (!strcmp("amd",options.ordering)) {
	     ierr=ZSYMperm_amd(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: amd performed\n");fflush(stdout);
#endif
	  }
	  else if (!strcmp("rcm",options.ordering)) {
	     ierr=ZSYMperm_rcm(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: rcm performed\n");fflush(stdout);
#endif
	  }
	  else { /* none */
	     ierr=ZSYMperm_null(A, rowscale,colscale, p,invq, &n, &options);
#ifdef PRINT_INFO
	     mexPrintf("ZSYMbilu: no symmetric reordering used\n");fflush(stdout);
#endif
	  }
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
	      pcolscale[i]=colscale[i].r;
	  } /* end for */
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif

	  /* no cosine blocking */
	  myblocksize=NULL;
	  mynblocks=0;
	  
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* set up the C for BILDL */
	  /* switch matrix A back to C notation and adapt data structures for bildl (and cosine_blocks) */
	  for (i=0; i<=n; i++)
	      A.ia[i]--;
	  for (i=0; i<A.ia[n]; i++)
	      A.ja[i]--;
	  /* use slightly different data structure for bildl */
	  C.nr=A.nr; C.nc=n;
	  C.ncol  =(integer *) malloc((size_t)n*sizeof(integer));
	  C.rowind=(integer **)malloc((size_t)n*sizeof(integer *));
	  C.val   =(doublecomplex **) malloc((size_t)n*sizeof(doublecomplex *));
	  C.nnz=A.nnz;
	  for (i=0; i<n; i++) {
	      C.ncol[i]  =A.ia[i+1]-A.ia[i];
	      C.rowind[i]=A.ja+A.ia[i];
	      C.val   [i]=A.a +A.ia[i];
	  } /* end for i */
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* we are ready to apply BILDL */
       }
    } /* end if no given blocks */
    else { /* block partitioning prescribed */
#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* do some scaling in advance */
	  ierr=ZSYMperm_null(A, rowscale,colscale, p,invq, &n, &options);
	  /* the scalings are real-valued anyway */
	  for (i=0; i<n; i++) {
	      pcolscale[i]=colscale[i].r;
	  } /* end for */
	  /* switch back to C-style */
	  for (i=0; i<n; i++) {
	      p[i]--; invq[i]--;
	  }
#ifdef _PROFILING_
	  time_reordering+=omp_get_wtime()-timeBegin;
#endif
#ifdef PRINT_CHECK
	  for (i=0;i<n; i++) {
	      if (p[i]<0 || p[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p[%ld]=%ld\n", i,i,p[i]);
		 fflush(stdout);
	      }
	      if (invq[i]<0 || invq[i]>=n) {
		 mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq[%ld]=%ld\n", i,i,invq[i]);
		 fflush(stdout);
	      }
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[p[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
		 fflush(stdout);
	      }
	      check_buff[p[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
	  for (i=0; i<n; i++) {
	      if (check_buff[invq[i]]) {
	         mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
		 fflush(stdout);
	      }
	      check_buff[invq[i]]=-1;
	  }
	  for (i=0; i<n; i++)
	      check_buff[i]=0;
#endif

	  myblocksize=options.blocksize;
	  mynblocks=options.nblocks;

#ifdef _PROFILING_
	  timeBegin=omp_get_wtime();
#endif
	  /* set up the C for BILDL */
	  /* switch matrix A back to C notation and adapt data structures for bildl (and cosine_blocks) */
	  for (i=0; i<=n; i++)
	      A.ia[i]--;
	  for (i=0; i<A.ia[n]; i++)
	      A.ja[i]--;
	  /* use slightly different data structure for bildl */
	  C.nr=A.nr; C.nc=n;
	  C.ncol  =(integer *)       malloc((size_t)n*sizeof(integer));
	  C.rowind=(integer **)      malloc((size_t)n*sizeof(integer *));
	  C.val   =(doublecomplex **)malloc((size_t)n*sizeof(doublecomplex *));
	  C.nnz=A.nnz;
	  for (i=0; i<n; i++) {
	      C.ncol[i]  =A.ia[i+1]-A.ia[i];
	      C.rowind[i]=A.ja+A.ia[i];
	      C.val   [i]=A.a +A.ia[i];
	  } /* end for i */
#ifdef _PROFILING_
	  time_pre_processing+=omp_get_wtime()-timeBegin;
#endif
	  /* before starting BILDL, we are going to reorder the system based on the
	     given block structure 
	  */
    } /* end if-else no given blocks */


    /* if only maximum weight matching was requested without cosine and if no
       initial block partitioning was passed (similar case for 
       no matching+no cosine+no blocks), then we are already done 
       otherwise (cosine or given blocks) we need to apply the symmetric 
       reordering to some compressed companion matrix
    */
    if (myblocksize!=NULL) {
#ifdef _PROFILING_
       timeBegin=omp_get_wtime();
#endif
       /* compute a companion matrix B of A(p,p) compressing the associated blocks */
       if (invblock==NULL)
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
       B.ia  =(integer *)malloc((size_t)(mynblocks+1)*sizeof(integer));
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
	       B.a[l].r=1.0; /* dummy */
	       B.a[l].i=0.0; /* dummy */
	       l++; 
	   } /* end for k */
	   cnt=0;  
       } /* end for i */  
       free(idx);
       free(idxpos);
#ifdef PRINT_INFO
       mexPrintf("ZSYMbilu: companion matrix for symmetric reorderings computed\n");fflush(stdout);
#endif

#ifdef PRINT_INFO2
       mexPrintf("ZSYMbilu: call (scaling+) reordering\n");fflush(stdout);
       for (i=0; i<B.nc; i++) {
	   mexPrintf("column %8ld\n",i+1);
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       k=B.ja[j-1];
	       mexPrintf("%8ld",k);
	   }
	   mexPrintf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       mexPrintf("%8.1le",B.a[j-1].r);
	   }
	   mexPrintf("\n");
	   for (j=B.ia[i]; j<B.ia[i+1]; j++) {
	       mexPrintf("%8.1le",B.a[j-1].i);
	   }
	   mexPrintf("\n");
       }
#endif

       /* select permutation and scaling driver by user request */
       nB=mynblocks;
       /* block permutation */
       p2   =(integer *)malloc((size_t)nB*sizeof(integer));
       invq2=(integer *)malloc((size_t)nB*sizeof(integer));
       
       if (!strcmp("mtmetis",options.ordering)) {
	  ierr=ZSYMperm_mtmetis(B, rowscale,colscale, p2,invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: mtmetis performed\n");fflush(stdout);
#endif
       }
#ifdef _USE_METIS4_
       else if (!strcmp("metisn",options.ordering)) {
	  ierr=ZSYMperm_metis_n(B, rowscale,colscale, p2,invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: metisn performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("metise",options.ordering)) {
	  ierr=ZSYMperm_metis_e(B, rowscale,colscale, p2,invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: metise performed\n");fflush(stdout);
#endif
       }
#endif // _USE_METIS4_
       else if (!strcmp("amd",options.ordering)) {
	  ierr=ZSYMperm_amd(B, rowscale,colscale, p2,invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: amd performed\n");fflush(stdout);
#endif
       }
       else if (!strcmp("rcm",options.ordering)) {
	  ierr=ZSYMperm_rcm(B, rowscale,colscale, p2,invq2, &nB, &options);
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: rcm performed\n");fflush(stdout);
#endif
       }
       else { /* none */
	  for (i=0; i<nB; i++)
	      p2[i]=invq2[i]=i+1;
#ifdef PRINT_INFO
	  mexPrintf("ZSYMbilu: no symmetric reordering used\n");fflush(stdout);
#endif
       }
      
    
#ifdef PRINT_INFO
       mexPrintf("ZSYMbilu: block (scaling +) reordering computed, ierr=%d\n",ierr);fflush(stdout);
#endif
#ifdef PRINT_INFO0
       mexPrintf("block p as returned\n");
       for (i=0; i<B.nc; i++) {
	   mexPrintf("%8ld",p2[i]);
       }
       mexPrintf("\n");
       fflush(stdout);
       mexPrintf("block invq as returned \n");
       fflush(stdout);
       for (i=0; i<B.nc; i++) {
	   mexPrintf("%8ld",invq2[i]);
       }
       mexPrintf("\n");
       fflush(stdout);
#endif
       free(B.ia);
       free(B.ja);
       free(B.a);
#ifdef PRINT_CHECK
       for (i=0;i<nB; i++) {
	   if (p2[i]<=0 || p2[i]>nB) {
	      mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, p2[%ld]=%ld\n", i,i,p2[i]);
	      fflush(stdout);
	   }
	   if (invq2[i]<=0 || invq2[i]>nB) {
	      mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq2[%ld]=%ld\n", i,i,invq2[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
       for (i=0; i<nB; i++) {
	   if (check_buff[p2[i]-1]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p2[%ld]=%ld\n",i,i,p2[i]);
	      fflush(stdout);
	   }
	   check_buff[p2[i]-1]=-1;
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
       for (i=0; i<nB; i++) {
	   if (check_buff[invq2[i]-1]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq2[%ld]=%ld\n",i,i,invq2[i]);
	      fflush(stdout);
	   }
	   check_buff[invq2[i]-1]=-1;
       }
       for (i=0; i<nB; i++)
	   check_buff[i]=0;
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
		 mexPrintf("ZSYMbilu ERROR, wrong block size, permuted cosine, step %ld, blocksize[%ld]=%ld\n", j,j,myblocksize[j]);
	      if (k>n)
		 mexPrintf("ZSYMbilu ERROR, cummulative block size exceeds system size, permuted cosine, step %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
	  }  
	  if (k!=n)
	     mexPrintf("ZSYMbilu ERROR, cummulative block does not match system size, permuted cosine, nblocks %ld, sum(blocksize)=%ld, n=%ld\n", j,k,n);
#ifdef PRINT_CHECK0
	  mexPrintf("myblocksize[0,...,%ld]\n",mynblocks-1);fflush(stdout);
	  for (j=0; j<mynblocks; j++)
	      mexPrintf("%4ld",myblocksize[j]);
	  mexPrintf("\n");fflush(stdout);
#endif
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
#ifdef PRINT_CHECK
       for (i=0;i<n; i++) {
	   if (invblock[i]<0 || invblock[i]>=n) {
	      mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invblock[%ld]=%ld\n", i,i,invblock[i]);
	      fflush(stdout);
	   }
	   if (invq3[i]<0 || invq3[i]>=n) {
	      mexPrintf("ZSYMbilu ERROR, wrong permutation product, step %ld, invq3[%ld]=%ld\n", i,i,invq3[i]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invblock[i]]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invblock[%ld]=%ld\n",i,i,invblock[i]);
	      fflush(stdout);
	   }
	   check_buff[invblock[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq3[i]]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq3[%ld]=%ld\n",i,i,invq3[i]);
	      fflush(stdout);
	   }
	   check_buff[invq3[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif

#ifdef PRINT_INFO
       mexPrintf("ZSYMbilu: permutation from reordering of the companion matrix expanded to full size\n");fflush(stdout);
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
#ifdef PRINT_CHECK0
       mexPrintf("p\n");
       for (j=0; j<n; j++)
	   mexPrintf("%6ld",p[j]);
       mexPrintf("\ninvq\n");
       fflush(stdout);
       for (j=0; j<n; j++)
	   mexPrintf("%6ld",invq[j]);
       mexPrintf("\n");
       fflush(stdout);
#endif
#ifdef PRINT_CHECK
       for (j=0; j<n; j++) {
	   if (invq[p[j]]!=j || p[j]<0 || p[j]>=n || invq[j]<0 || invq[j]>=n) {
	      mexPrintf("ZSYMbilu ERROR, wrong p/invq after reorering the compressed graph, step %ld, p[%ld]=%ld, invq[%ld]=%ld\n", j,j,p[j],p[j],invq[p[j]]);
	      fflush(stdout);
	   }
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[p[i]]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, p[%ld]=%ld\n",i,i,p[i]);
	      fflush(stdout);
	   }
	   check_buff[p[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
       for (i=0; i<n; i++) {
	   if (check_buff[invq[i]]) {
	      mexPrintf("ZSYMbilu ERROR, duplicate index, step %ld, invq[%ld]=%ld\n",i,i,invq[i]);
	      fflush(stdout);
	   }
	   check_buff[invq[i]]=-1;
       }
       for (i=0; i<n; i++)
	   check_buff[i]=0;
#endif

#ifdef PRINT_INFO
       mexPrintf("ZSYMbilu: permutation product built\n");fflush(stdout);
#endif
    
       free(invblock);
       free(p3);
       free(invq3);
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
       mexPrintf("ZSYMbilu: matrix for BILU prepared\n");fflush(stdout);
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
    } /* end if blocksize!=0 */
#ifdef _PROFILING_
    time_reordering+=omp_get_wtime()-timeBegin;
#endif
    

      

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
    /* start already exporting the scalings */
    /* 4. lhs: S_L */
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
    mexPrintf("ZSYMbilu: S_L already exported\n");fflush(stdout);
#endif

    /* pcolscale not required anymore since C is already properly scaled */
    free(colscale);
    free(pcolscale);
#ifdef _PROFILING_
    time_post_processing=omp_get_wtime()-timeBegin;
#endif
    
#ifdef PRINT_INFO
    mexPrintf("call SYMBILU, droptol=%8.1le\n", options.droptol); fflush(stdout);
#endif

    

#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif   
    /* blocksize is now incorporated by the permutations */
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
    /* mexPrintf("ilu1t=%ld, mynblocks=%ld, blocksize=%ld\n",ilu1t, mynblocks,myblocksize); */
    if (!invert_blocks)
       pivots=(integer *)malloc((size_t)n*sizeof(integer));
    else
       pivots=NULL;
    ierr=ZSYMbilu(&C, &BL,&BiD,&BL, NULL,NULL, p,invq, &determinant, &isdefinite,
		  myblocksize,mynblocks,
		  options.droptol,ilu1t,pa,perturbation,pivots,level_of_fill);
    for (i=0; i<n; i++)
        determinant.r-=2*log(valuesR[i]);

#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: bilu completed, ierr=%d\n",ierr);fflush(stdout);
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
       mexErrMsgTxt("ZSYMbilu: block ILU abnormally terminated") ;
    }
#ifdef _PROFILING_
    time_bilu=omp_get_wtime()-timeBegin;
#endif


#ifdef _PROFILING_
    timeBegin=omp_get_wtime();
#endif
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
    free(invq);
    free(p);
#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: P exported\n");fflush(stdout);
#endif

#ifdef PRINT_INFO
    if (options.blocksize==NULL)
       mexPrintf("ZSYMbilu: options.blocksize unused\n");
    else
       mexPrintf("ZSYMbilu: export options.blocksize\n");
    fflush(stdout);
#endif
    
    /* 4. lhs: determinant */
    plhs[4]=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxCOMPLEX);
    pr=mxGetPr(plhs[4]);
    *pr=determinant.r;
    pr=mxGetPi(plhs[4]);
    *pr=determinant.i;
    
    /* 5. lhs: pivots, if used */
    if (invert_blocks)
       plhs[5]=mxCreateDoubleMatrix((mwSize)1,(mwSize)0, mxREAL);
    else {
       plhs[5]=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
       pr=mxGetPr(plhs[5]);
       for (i=0; i<n; i++)
	   pr[i]=pivots[i];
       free(pivots);
    } 
    
    /* export options */
    if (options.blocksize!=NULL) {
       plhs[6]=mxCreateStructMatrix((mwSize)1, (mwSize)1, 11, optionsnamesx);
    }
    else {
       plhs[6]=mxCreateStructMatrix((mwSize)1, (mwSize)1, 10, optionsnames);
    }
    options_output=plhs[6];
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
    /* field 3: matching */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=options.matching;
    mxSetFieldByNumber(options_output, (mwSize)0, 3, tmp);
    /* field 5: droptol */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=options.droptol;
    mxSetFieldByNumber(options_output, (mwSize)0, 5, tmp);
    /* field 4: ordering */    
    output_buf=(char *)mxCalloc((size_t)strlen(options.ordering)+1, (size_t)sizeof(char));
    strcpy(output_buf,options.ordering);
    if (flag)
       free(options.ordering);
    fout = mxCreateString(output_buf);
    mxSetFieldByNumber(options_output, (mwSize)0, 4, fout);
    mxFree(output_buf);
    /* field 6: perturbation */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=perturbation;
    mxSetFieldByNumber(options_output, (mwSize)0, 6, tmp);
    /* field 7: isdefinite */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=isdefinite;
    mxSetFieldByNumber(options_output, (mwSize)0, 7, tmp);
    /* field 8: invert_blocks */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=invert_blocks;
    mxSetFieldByNumber(options_output, (mwSize)0, 8, tmp);
    /* field 9: level_of_fill */
    tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(tmp);
    *pr=level_of_fill;
    mxSetFieldByNumber(options_output, (mwSize)0, 9, tmp);
    /* field 10: blocksize (optionally) */    
    if (options.blocksize!=NULL) {
       tmp=mxCreateDoubleMatrix((mwSize)1,(mwSize)options.nblocks, mxREAL);
       pr=mxGetPr(tmp);
       for (i=0; i<options.nblocks; i++) {
	   pr[i]=options.blocksize[i];
       } /* end for i */
       free(options.blocksize);
       mxSetFieldByNumber(options_output, (mwSize)0, 10, tmp);
    }
    
#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: options structure exported\n");fflush(stdout);
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
    mexPrintf("ZSYMbilu: BL exported\n");fflush(stdout);
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
	   /* scalar case: pass psychological "D factor" dii rather than 1/dii */
	   if (m==1) {
	      val=sqrt(pr[0]*pr[0]+pi[0]*pi[0]);
	      pr[0]= pr[0]/val;
	      pi[0]=-pi[0]/val;
	   }
	   else { /* set strict upper triangular part to 0 */
	      pr+=m;
	      for (j=1; j<m; j++) {
		  for (k=0; k<j; k++)
		      pr[k]=pi[k]=0.0;
		  pr+=m;
	      } /* end for k */
	   } /* end if-else m=1 */
	} /* end if */
	/* set each field in Block structure */
	mxSetFieldByNumber(Block, (mwIndex)0, 3, BlockD);

	/* finally set output BiD{i} */
	mxSetCell(BiD_output, (mwIndex)i,Block);
    } /* end for i */
#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: BiD exported\n");fflush(stdout);
#endif

    /* finally release auxiliary memory */
    ZSparseBlockDelete(&BL);
    ZSparseBlockDelete(&BiD);
#ifdef _PROFILING_
    time_post_processing+=omp_get_wtime()-timeBegin;
    time_total=omp_get_wtime()-time_total;
    mexPrintf("CMEX profiling summary\n");
    mexPrintf("  i) pre processing         %8.1le\n",time_pre_processing);
    mexPrintf(" ii) cosine-based blocking  %8.1le\n",time_cosine);
    mexPrintf("iii) reordering             %8.1le\n",time_reordering);
    mexPrintf(" iv) BILDL driver           %8.1le\n",time_bilu);
    mexPrintf("  v) post processing        %8.1le\n",time_post_processing);
    mexPrintf("Total CMEX-BILDL time %8.1le\n\n",time_total);

    fflush(stdout);
#endif


#ifdef PRINT_INFO
    mexPrintf("ZSYMbilu: memory released\n");fflush(stdout);
#endif
#ifdef PRINT_CHECK
    free(check_buff);
#endif
    return;       

}
