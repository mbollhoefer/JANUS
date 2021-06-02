/* $Id: janussolver.c 6094 2020-03-28 16:17:50Z bolle $

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	August 29, 2017. JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define STDERR          stderr
#define STDOUT          stdout


// #define PRINT_INFO
// #define PRINT_CHECK
#define _PROFILING_

#ifdef _PROFILING_
#include <omp.h>
#endif
integer JanusSolver(SparseMatrix *A, JanusPrec *PREC, void *rhs, void *sol,
		    integer K_RESTART, double RES_TOL, integer MAX_IT, integer *iter)
{
   integer ierr=0, ipar[16],flag, n=PREC->n, i,j,k;
   double  fpar[16], v, nrm1, nrm, nrm1inf, nrminf;


#ifdef PRINT_INFO
   printf("entering iterative solver\n");fflush(stdout);
#endif
   
#ifdef _PROFILING_
   double timeBegin,
          time_krylov=0.0;
#endif
  
#ifdef _PROFILING_
   timeBegin=omp_get_wtime();
#endif
   if (A->isreal) {
#ifdef PRINT_INFO
      printf("entering real-valued solver\n");fflush(stdout);
#endif
      double  *drhs=(double *) rhs, *dsol=(double *)sol, *dw, *dy,
	      **pdv=(double **)A->val;
      dy=(double *) malloc((size_t)n*sizeof(double));
      
      if (A->issymmetric|| A->ishermitian) {
	 // SPD case: preconditioned CG
	 if (A->isdefinite) {
#ifdef PRINT_INFO
	    printf("entering real-symmetric positive definite solver\n");fflush(stdout);
#endif
	    dw=(double *) malloc((size_t)5*n*sizeof(double));
	    // CG set up 
	    ipar[0]=0;      // first call of cg
	    ipar[1]=2;      // use preconditioning
	    ipar[2]=3;      // relative error in the energy norm as stopping criterion
	    ipar[3]=n*5;    // size of w
	    ipar[5]=MAX_IT; // maximum number of iteration steps
	    ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

	    fpar[0]=RES_TOL;// relative error tolerance
	    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[11]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

	    // main loop using reverse communication
	    flag=-1;
	    while (flag) {
	          Dpcg(&n,drhs,dsol,ipar,fpar,dw);

		  // what is required by the solver?
		  dw--;
		  switch (ipar[0]) {
		  case  1:   // apply matrix vector multiplication
		  case  2:
		        MatVec(A,(void *)(dw+ipar[7]),(void *)(dw+ipar[8]));
			break;
		  case  3:   // apply preconditioner
		  case  4:   
		  case  5:   
		  case  6:   
		        JanusSol(PREC, (void *)(dw+ipar[7]),(void *)(dw+ipar[8]),dy,1);
			break;
		  default:   // the iterative solver terminated with code=ipar[0]
		        ierr=ipar[0];
			flag=0;
			break;
		  } // end switch
		  dw++;
	    } // end while

	    
	    // why did CG stop?
	    // error?
	    if (ierr) {
	       if (ierr==-1)
		  printf("too many iteration steps\n");
	       else if (ierr==-2)
		  printf("not enough work space\n");
	       else if (ierr==-3)
		  printf("algorithm breaks down\n");
	       else 
		  printf("algorithm stops with error code %ld\n",ierr);
	    }
	    else { // success
	       if (ipar[6]>=MAX_IT)
		  printf("preconditioned CG stopped after %ld iteration steps\n", ipar[6]);
	       else
		  printf("preconditioned CG successfully completed with %ld iteration steps\n", ipar[6]);
	    }
	    // ------------------- END CG preconditioned by BICHOL --------------------
	    // ------------------------------------------------------------------------
	 }
	 // symmetric indefinite case: preconditioned SQMR
	 else { // non-SPD case
#ifdef PRINT_INFO
	    printf("entering real-symmetric solver\n");fflush(stdout);
#endif
	    dw=(double *) malloc((size_t)6*n*sizeof(double));
	    // SQMR set up

	    // since I prefer the backward error as stopping criterion, I provide
	    // an estimate for ||A||_2
	    // use ||A||_inf=||A||_1 as estimate
	    // use w as buffer for the matrix 1/inf norm
	    for (i=0; i<n; i++)
	       dw[i]=0.0;
	    for (j=0; j<n; j++) {
	        for (k=0; k<A->ncol[j]; k++) {
		    // row index i of A(i,j)
		    i=A->rowind[j][k];
		    // numerical value A(i,j)
		    v=fabs(pdv[j][k]);
		    // contribute to 1-norm of A(:,i)
		    dw[i]+=v;
		    // contribute to 1-norm of A(:,j) (half of the matrix is missing)
		    if (i!=j)
		       dw[j]+=v;
		} // end for k
	    } // end for i
	    // build matrix 1/inf-norm
	    nrm1inf=0.0;
	    for (i=0; i<n; i++) {
	        nrm1inf=MAX(nrm1inf,dw[i]);
	    } // end for i
	    // printf("%8.1le\n",nrm1inf);
	    // end computation of ||A||_1=||A||_inf
    
	    ipar[0]=0;       // first call of sqmr
	    ipar[1]=2;       // use preconditioning
	    ipar[2]=3;       // backward error as stopping criterion
	    ipar[3]=n*6;     // size of w
	    ipar[5]=MAX_IT;  // maximum number of iteration steps
	    ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

	    fpar[0]=RES_TOL; // relative error tolerance
	    fpar[11]=nrm1inf;// estimate for the matrix-2 norm
	    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

	    // main loop using reverse communication
	    flag=-1;
	    while (flag) {
	          DSYMqmr(&n,drhs,dsol,ipar,fpar,dw);

		  // what is required by the solver?
		  dw--;
		  // printf("ipar[0]=%ld, ipar[6]=%ld, fpar[3]=%8.1le, fpar[4]=%8.1le\n",ipar[0],ipar[6],fpar[3],fpar[4]);
		  switch (ipar[0]) {
		  case  1:   // apply matrix vector multiplication
		  case  2:
		        MatVec(A,(void *)(dw+ipar[7]),(void *)(dw+ipar[8]));
			break;
		  case  3:   // apply preconditioner
		  case  4:   
		  case  5:   
		  case  6:   
		        JanusSol(PREC, (void *)(dw+ipar[7]),(void *)(dw+ipar[8]),dy,1);
			break;
		  default:   // the iterative solver terminated with code=ipar[0]
		    ierr=ipar[0];
		    flag=0;
		    break;
		  } // end switch
		  dw++;
	    } // end while

	    // why did SQMR stop?
	    // error?
	    if (ierr) {
	       if (ierr==-1)
		  printf("too many iteration steps\n");
	       else if (ierr==-2)
		  printf("not enough work space\n");
	       else if (ierr==-3)
		  printf("algorithm breaks down\n");
	       else 
		  printf("algorithm stops with error code %ld\n",ierr);
	    }
	    else { // success
	       if (ipar[6]>=MAX_IT)
		  printf("preconditioned SQMR stopped after %ld iteration steps\n", ipar[6]);
	       else
		  printf("preconditioned SQMR successfully completed with %ld iteration steps\n", ipar[6]);
	    }
	    // ------------------- END SQMR preconditioned by BILDL ------------------
	    // ------------------------------------------------------------------------
	 }
      }
      else { // non-symmetric case
#ifdef PRINT_INFO
	 printf("entering real general solver\n");fflush(stdout);
#endif
	 dw=(double *)malloc((size_t)((n+3)*(K_RESTART+2)+((K_RESTART+1)*K_RESTART)/2)*sizeof(double));

	 // GMRES set up

	 // since I prefer the backward error as stopping criterion, I provide
	 // an estimate for ||A||_2
	 // use geometric mean of ||A||_inf and ||A||_1 as estimate
	 // use w as buffer for the matrix inf norm
	 for (i=0; i<n; i++)
	     dw[i]=0.0;
	 nrm1=0.0;
	 for (j=0; j<n; j++) {
	     nrm=0.0;
	     for (k=0; k<A->ncol[j]; k++) {
	         // row index i of A(i,j)
	         i=A->rowind[j][k];
		 // numerical value A(i,j)
		 v=fabs(pdv[j][k]);
		 // compute 1-norm of A(:,j)
		 nrm+=v;
		 // contribute to 1-norm of A(i,:)
		 dw[i]+=v;
	     } // end for k
	     // build matrix 1-norm
	     nrm1=MAX(nrm1,nrm);
	 } // end for i
	 nrminf=0.0;
	 for (i=0; i<n; i++) {
	     // build matrix inf-norm
	     nrminf=MAX(nrminf,dw[i]);
	 } // end for i
	 // end computation of ||A||_1, ||A||_inf

    
	 ipar[0]=0;        // first call of GMRES
	 ipar[1]=2;        // use right preconditioning
	 ipar[2]=3;        // relative error in the energy norm as stopping criterion
	 ipar[3]=(n+3)*(K_RESTART+2)+((K_RESTART+1)*K_RESTART)/2;    // size of w
	 ipar[4]=K_RESTART;// restart length
	 ipar[5]=MAX_IT;   // maximum number of iteration steps
	 ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

	 fpar[0]=RES_TOL; // relative error tolerance
	 fpar[11]=sqrt(nrm1*nrminf);// estimate for the matrix-2 norm
	 // printf("A-nrm=%8.1le\n",fpar[11]);
	 fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

	 // main loop using reverse communication
	 flag=-1;
	 while (flag) {
	       Dgmres(&n,drhs,dsol,ipar,fpar,dw);

	       // printf("%ld\n",ipar[0]);
	       // what is required by the solver?
	       dw--;
	       switch (ipar[0]) {
	       case  1:   // apply matrix vector multiplication
	       case  2:
		     // printf("%4ld, %8.1le, %8.1le\n",ipar[6],fpar[4],fpar[3]);
		     MatVec(A,(void *)(dw+ipar[7]),(void *)(dw+ipar[8]));
		     break;
	       case  3:   // apply preconditioner
	       case  4:   
	       case  5:   
	       case  6:   
		     JanusSol(PREC, (void *)(dw+ipar[7]),(void *)(dw+ipar[8]),dy,1);
		     break;
	       default:   // the iterative solver terminated with code=ipar[0]
		     ierr=ipar[0];
		     flag=0;
		     break;
	       } // end switch
	       dw++;
	 } // end while

	 // why did GMRES stop?
	 // error?
	 if (ierr) {
	    if (ierr==-1)
	       printf("too many iteration steps\n");
	    else if (i==-2)
	       printf("not enough work space\n");
	    else if (i==-3)
	       printf("algorithm breaks down\n");
	    else 
	       printf("algorithm stops with error code %ld\n",ierr);
	 }
	 else { // success
	    if (ipar[6]>=MAX_IT)
	       printf("preconditioned GMRES stopped after %ld iteration steps\n", ipar[6]);
	    else
	       printf("preconditioned GMRES successfully completed with %ld iteration steps\n", ipar[6]);
	 }
      }

      free(dy);
      free(dw);
   }
   else { //non-real case
#ifdef PRINT_INFO
      printf("entering complex-valued solver\n");fflush(stdout);
#endif
      doublecomplex  *zrhs=(doublecomplex *) rhs, *zsol=(doublecomplex *)sol, *zw, *zy,
	             **pzv=(doublecomplex **)A->val;
      zy=(doublecomplex *) malloc((size_t)n*sizeof(doublecomplex));

      if (A->issymmetric) {
#ifdef PRINT_INFO
	 printf("entering complex-symmetric solver\n");fflush(stdout);
#endif
	 zw=(doublecomplex *) malloc((size_t)6*n*sizeof(doublecomplex));
	 
	 // SQMR set up 

	 // since I prefer the backward error as stopping criterion, I provide
	 // an estimate for ||A||_2
	 // use ||A||_inf=||A||_1 as estimate
	 // use w as buffer for the matrix 1/inf norm
	 for (i=0; i<n; i++)
	     zw[i].r=0.0;
	 for (j=0; j<n; j++) {
	     for (k=0; k<A->ncol[j]; k++) {
	         // row index i of A(i,j)
	         i=A->rowind[j][k];
		 // numerical value A(i,j)
		 v=sqrt(pzv[j][k].r*pzv[j][k].r+pzv[j][k].i*pzv[j][k].i);
		 // contribute to 1-norm of A(:,i)
		 zw[i].r+=v;
		 // contribute to 1-norm of A(:,j) (half of the matrix is missing)
		 if (i!=j)
		    zw[j].r+=v;
	     } // end for k
	 } // end for i
	 // build matrix 1/inf-norm
	 nrm1inf=0.0;
	 for (i=0; i<n; i++) {
	     nrm1inf=MAX(nrm1inf,zw[i].r);
	 } // end for i
	 // printf("%8.1le\n",nrm1inf);
	 // end computation of ||A||_1=||A||_inf
    
	 ipar[0]=0;      // first call of cg
	 ipar[1]=2;      // use preconditioning
	 ipar[2]=3;      // relative error in the energy norm as stopping criterion
	 ipar[3]=n*6;    // size of w
	 ipar[5]=MAX_IT; // maximum number of iteration steps
	 ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;
	 
	 fpar[0]=RES_TOL;// relative error tolerance
	 fpar[11]=nrm1inf;// estimate for the matrix-2 norm
	 fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;
	 
	 // main loop using reverse communication
	 flag=-1;
	 while (flag) {
#ifdef PRINT_INFO
	       printf("call SQMR\n");fflush(stdout);
#endif
	       ZSYMqmr(&n,zrhs,zsol,ipar,fpar,zw);

	       // what is required by the solver?
	       zw--;
	       // printf("ipar[0]=%ld, ipar[6]=%ld, fpar[3]=%8.1le, fpar[4]=%8.1le\n",ipar[0],ipar[6],fpar[3],fpar[4]);
	       switch (ipar[0]) {
	       case  1:   // apply matrix vector multiplication
	       case  2:
#ifdef PRINT_INFO
		     printf("MatVec\n");fflush(stdout);
#endif
		     MatVec(A,(void *)(zw+ipar[7]),(void *)(zw+ipar[8]));
		     break;
	       case  3:   // apply preconditioner
	       case  4:   
	       case  5:   
	       case  6:   
#ifdef PRINT_INFO
		     printf("Forward/Backward Solve\n");fflush(stdout);
#endif
		     JanusSol(PREC, (void *)(zw+ipar[7]),(void *)(zw+ipar[8]),zy,1);
		     break;
	       default:   // the iterative solver terminated with code=ipar[0]
		     ierr=ipar[0];
		     flag=0;
		     break;
	       } // end switch
	       zw++;
	 } // end while

	 // why did SQMR stop?
	 // error?
	 if (ierr) {
	    if (ierr==-1)
	       printf("too many iteration steps\n");
	    else if (ierr==-2)
	       printf("not enough work space\n");
	    else if (ierr==-3)
	       printf("algorithm breaks down\n");
	    else 
	       printf("algorithm stops with error code %ld\n",ierr);
	 }
	 else { // success
 	    if (ipar[6]>=MAX_IT)
	       printf("preconditioned SQMR stopped after %ld iteration steps\n", ipar[6]);
	    else
 	       printf("preconditioned SQMR successfully completed with %ld iteration steps\n", ipar[6]);
	 }
	 // ------------------- END SQMR preconditioned by BILDL --------------------
	 // ------------------------------------------------------------------------
      }
      else if (A->ishermitian) {
#ifdef PRINT_INFO
	 printf("entering complex-Hermitian positive definite solver\n");fflush(stdout);
#endif
	 // HPD case: preconditioned CG
	 if (A->isdefinite) {
	    zw=(doublecomplex *) malloc((size_t)5*n*sizeof(doublecomplex));
	    // CG set up 
	    ipar[0]=0;      // first call of cg
	    ipar[1]=2;      // use preconditioning
	    ipar[2]=3;      // relative error in the energy norm as stopping criterion
	    ipar[3]=n*5;    // size of w
	    ipar[5]=MAX_IT; // maximum number of iteration steps
	    ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

	    fpar[0]=RES_TOL;// relative error tolerance
	    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[11]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

	    // main loop using reverse communication
	    flag=-1;
	    while (flag) {
#ifdef PRINT_INFO
	          printf("Call CG\n");fflush(stdout);
#endif
	          Zpcg(&n,zrhs,zsol,ipar,fpar,zw);

		  // what is required by the solver?
		  zw--;
		  switch (ipar[0]) {
		  case  1:   // apply matrix vector multiplication
		  case  2:
#ifdef PRINT_INFO
		        printf("MatVec\n");fflush(stdout);
#endif
		        MatVec(A,(void *)(zw+ipar[7]),(void *)(zw+ipar[8]));
			break;
		  case  3:   // apply preconditioner
		  case  4:   
		  case  5:   
		  case  6:   
#ifdef PRINT_INFO
		        printf("Forward/Backward Solve\n");fflush(stdout);
#endif
		        JanusSol(PREC, (void *)(zw+ipar[7]),(void *)(zw+ipar[8]),zy,1);
			break;
		  default:   // the iterative solver terminated with code=ipar[0]
		        ierr=ipar[0];
			flag=0;
			break;
		  } // end switch
		  zw++;
	    } // end while

	    
	    // why did CG stop?
	    // error?
	    if (ierr) {
	       if (ierr==-1)
		  printf("too many iteration steps\n");
	       else if (ierr==-2)
		  printf("not enough work space\n");
	       else if (ierr==-3)
		  printf("algorithm breaks down\n");
	       else 
		  printf("algorithm stops with error code %ld\n",ierr);
	    }
	    else { // success
	       if (ipar[6]>=MAX_IT)
		  printf("preconditioned CG stopped after %ld iteration steps\n", ipar[6]);
	       else
		  printf("preconditioned CG successfully completed with %ld iteration steps\n", ipar[6]);
	    }
	 }
	 else {
#ifdef PRINT_INFO
	    printf("entering complex-Hermitian solver\n");fflush(stdout);
#endif
	    zw=(doublecomplex *) malloc((size_t)6*n*sizeof(doublecomplex));
	 
	    // SQMR set up 

	    // since I prefer the backward error as stopping criterion, I provide
	    // an estimate for ||A||_2
	    // use ||A||_inf=||A||_1 as estimate
	    // use w as buffer for the matrix 1/inf norm
	    for (i=0; i<n; i++)
	        zw[i].r=0.0;
	    for (j=0; j<n; j++) {
	        for (k=0; k<A->ncol[j]; k++) {
		    // row index i of A(i,j)
		    i=A->rowind[j][k];
		    // numerical value A(i,j)
		    v=sqrt(pzv[j][k].r*pzv[j][k].r+pzv[j][k].i*pzv[j][k].i);
		    // contribute to 1-norm of A(:,i)
		    zw[i].r+=v;
		    // contribute to 1-norm of A(:,j) (half of the matrix is missing)
		    if (i!=j)
		       zw[j].r+=v;
		} // end for k
	    } // end for i
	    // build matrix 1/inf-norm
	    nrm1inf=0.0;
	    for (i=0; i<n; i++) {
	        nrm1inf=MAX(nrm1inf,zw[i].r);
	    } // end for i
	    // printf("%8.1le\n",nrm1inf);
	    // end computation of ||A||_1=||A||_inf
    
	    ipar[0]=0;      // first call of sqmr
	    ipar[1]=2;      // use preconditioning
	    ipar[2]=3;      // relative error in the energy norm as stopping criterion
	    ipar[3]=n*6;    // size of w
	    ipar[5]=MAX_IT; // maximum number of iteration steps
	    ipar[4]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;
	 
	    fpar[0]=RES_TOL;// relative error tolerance
	    fpar[11]=nrm1inf;// estimate for the matrix-2 norm
	    fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;
	 
	    // main loop using reverse communication
	    flag=-1;
	    while (flag) {
#ifdef PRINT_INFO
	          printf("Call HQMR\n");fflush(stdout);
#endif
	          ZHERqmr(&n,zrhs,zsol,ipar,fpar,zw);

		  // what is required by the solver?
		  zw--;
		  // printf("ipar[0]=%ld, ipar[6]=%ld, fpar[3]=%8.1le, fpar[4]=%8.1le\n",ipar[0],ipar[6],fpar[3],fpar[4]);
		  switch (ipar[0]) {
		  case  1:   // apply matrix vector multiplication
		  case  2:
#ifdef PRINT_INFO
		        printf("MatVec\n");fflush(stdout);
#endif
		        MatVec(A,(void *)(zw+ipar[7]),(void *)(zw+ipar[8]));
			break;
		  case  3:   // apply preconditioner
		  case  4:   
		  case  5:   
		  case  6:   
#ifdef PRINT_INFO
		        printf("Forward/Backward Solve\n");fflush(stdout);
#endif
		        JanusSol(PREC, (void *)(zw+ipar[7]),(void *)(zw+ipar[8]),zy,1);
			break;
		  default:   // the iterative solver terminated with code=ipar[0]
		        ierr=ipar[0];
			flag=0;
			break;
		  } // end switch
		  zw++;
	    } // end while

	    // why did SQMR stop?
	    // error?
	    if (ierr) {
	       if (ierr==-1)
		  printf("too many iteration steps\n");
	       else if (ierr==-2)
		  printf("not enough work space\n");
	       else if (ierr==-3)
		  printf("algorithm breaks down\n");
	       else 
		  printf("algorithm stops with error code %ld\n",ierr);
	    }
	    else { // success
	       if (ipar[6]>=MAX_IT)
		  printf("preconditioned SQMR stopped after %ld iteration steps\n", ipar[6]);
	       else
		  printf("preconditioned SQMR successfully completed with %ld iteration steps\n", ipar[6]);
	    }
	    // ------------------- END SQMR preconditioned by BILDL --------------------
	    // ------------------------------------------------------------------------
	 }
      }
      else { // general non-real case
#ifdef PRINT_INFO
	 printf("entering complex general solver\n");fflush(stdout);
#endif
	 zw=(doublecomplex *)malloc((size_t)((n+3)*(K_RESTART+2)+((K_RESTART+1)*K_RESTART)/2)*sizeof(doublecomplex));

	 // GMRES set up

	 // since I prefer the backward error as stopping criterion, I provide
	 // an estimate for ||A||_2
	 // use geometric mean of ||A||_inf and ||A||_1 as estimate
	 // use w as buffer for the matrix inf norm
	 for (i=0; i<n; i++)
	     zw[i].r=0.0;
	 nrm1=0.0;
	 for (j=0; j<n; j++) {
	     nrm=0.0;
	     for (k=0; k<A->ncol[j]; k++) {
	         // row index i of A(i,j)
	         i=A->rowind[j][k];
		 // numerical value A(i,j)
		 v=sqrt(pzv[j][k].r*pzv[j][k].r+pzv[j][k].i*pzv[j][k].i);
		 // compute 1-norm of A(:,j)
		 nrm+=v;
		 // contribute to 1-norm of A(i,:)
		 zw[i].r+=v;
	     } // end for k
	     // build matrix 1-norm
	     nrm1=MAX(nrm1,nrm);
	 } // end for i
	 nrminf=0.0;
	 for (i=0; i<n; i++) {
	     // build matrix inf-norm
	     nrminf=MAX(nrminf,zw[i].r);
	 } // end for i
	 // end computation of ||A||_1, ||A||_inf

    
	 ipar[0]=0;        // first call of GMRES
	 ipar[1]=2;        // use right preconditioning
	 ipar[2]=3;        // relative error in the energy norm as stopping criterion
	 ipar[3]=(n+3)*(K_RESTART+2)+((K_RESTART+1)*K_RESTART)/2;    // size of w
	 ipar[4]=K_RESTART;// restart length
	 ipar[5]=MAX_IT;   // maximum number of iteration steps
	 ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=ipar[14]=ipar[15]=0;

	 fpar[0]=RES_TOL; // relative error tolerance
	 fpar[11]=sqrt(nrm1*nrminf);// estimate for the matrix-2 norm
	 fpar[1]=fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

	 // main loop using reverse communication
	 flag=-1;
	 while (flag) {
#ifdef PRINT_INFO
	       printf("Call GMRES\n");fflush(stdout);
#endif
	       Zgmres(&n,zrhs,zsol,ipar,fpar,zw);

	       // printf("%ld\n",ipar[0]);
	       // what is required by the solver?
	       zw--;
	       switch (ipar[0]) {
	       case  1:   // apply matrix vector multiplication
	       case  2:
#ifdef PRINT_INFO
		     printf("MatVec\n");fflush(stdout);
#endif
		     // printf("%4ld, %8.1le, %8.1le\n",ipar[6],fpar[4],fpar[3]);
		     MatVec(A,(void *)(zw+ipar[7]),(void *)(zw+ipar[8]));
		     break;
	       case  3:   // apply preconditioner
	       case  4:   
	       case  5:   
	       case  6:   
#ifdef PRINT_INFO
		     printf("Forward/Backward Solve\n");fflush(stdout);
#endif
		     JanusSol(PREC, (void *)(zw+ipar[7]),(void *)(zw+ipar[8]),zy,1);
		     break;
	       default:   // the iterative solver terminated with code=ipar[0]
		     ierr=ipar[0];
		     flag=0;
		     break;
	       } // end switch
	       zw++;
	 } // end while

	 // why did GMRES stop?
	 // error?
	 if (ierr) {
	    if (ierr==-1)
	       printf("too many iteration steps\n");
	    else if (i==-2)
	       printf("not enough work space\n");
	    else if (i==-3)
	       printf("algorithm breaks down\n");
	    else 
	       printf("algorithm stops with error code %ld\n",ierr);
	 }
	 else { // success
	    if (ipar[6]>=MAX_IT)
	       printf("preconditioned GMRES stopped after %ld iteration steps\n", ipar[6]);
	    else
	       printf("preconditioned GMRES successfully completed with %ld iteration steps\n", ipar[6]);
	 }
      }

      free(zy);
      free(zw);
   }

   
#ifdef _PROFILING_
    time_krylov=omp_get_wtime()-timeBegin;
    // printf("Krylov subspace solver %8.1le[sec]\n",time_krylov);fflush(stdout);
#endif


   *iter=ipar[6];
   return (ierr);
}
