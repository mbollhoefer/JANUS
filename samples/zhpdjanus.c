#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <janus.h>


int main(int argc, char **argv)
{
    // prescribe lower triangular part of a sparse Hermitian pos. def. matrix
    long int      ncol[8]={  4,    3,    2,    2,   3,   2,   1,   1},
                  rowind0[4]={0,2,5,6},
	          rowind1[3]=     {1,2,4},
		  rowind2[2]=           {2,7},
		  rowind3[2]=                 {3,6},
		  rowind4[3]=                      {4,5,6},
		  rowind5[2]=                           {5,7},
		  rowind6[1]=                                {6},
		  rowind7[1]=                                     {7},
		  *rowind[8]={rowind0,rowind1,rowind2,rowind3,rowind4,rowind5,rowind6,rowind7};
    doublecomplex val0[4]   ={{17.0,0.0},{0.5,0.5},{2.0,0.0},{0.0,5.0}},
                  val1[3]   ={{8.0,0.0},{-3.0,1.0},{0.0,-2.0}},
		  val2[2]   ={{7.0,0.0},{1.0,0.0}},
		  val3[2]   ={{3.0,0.0},{0.0,1.0}},
		  val4[3]   ={{6.0,0.0},{0.5,0.5},{0.0,-5.0}},
		  val5[2]   ={{6.0,0.0},{4.0,1.0}},
		  val6[1]   ={{11.0,0.0}},
		  val7[1]   ={{7.0,0.0}},
                  *val[8]   ={val0,val1,val2,val3,val4,val5,val6,val7};

    long int nblocks,j,mx,iter,ierr;
    double   fill_bilu,mu,stddev;

    doublecomplex sol[8]={ {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}, {0.0,0.0}},
                  rhs[8]={ {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}, {1.0,0.0}};

    SparseMatrix  A;
    JanusOptions options;
    JanusPrec PREC;

    
    /* initialize matrix */
    A.nr=A.nc=8;
    A.nnz=18;
    A.ncol=ncol;
    A.rowind=rowind;
    A.val=(void *)val;
    /* tag properties of A */
    A.isreal=0;
    A.issymmetric=0;
    A.ishermitian=1;
    A.isdefinite=1;
    
    PrintMatrix(&A);



    //-------------------------------------------------------------------------
    //--------------- set up default user options for block ICHOL -------------
    //-------------------------------------------------------------------------
    JanusDefaultOptions(&options);



    //-------------------------------------------------------------------------
    //---------------------------- compute block ICHOL ------------------------
    //-------------------------------------------------------------------------
    ierr=JanusFactor(&A, &PREC, options);
    if (ierr) {
       printf("computation of block incomplete CHOL failed\n");
       JanusDelete(&PREC);
       return (-1);
    }


    
    //-------------------------------------------------------------------------
    // ------ Krylov subspace solver (here: CG) preconditioned by BICHOL ------
    // ------------------------------------------------------------------------
    ierr=JanusSolver(&A,&PREC,rhs,sol,30,1e-6,5000,&iter);
    // why did Krylov subspace solver stop?
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
       if (iter>=5000)
	  printf("preconditioned Krylov subspace solver stopped after %ld iteration steps\n", iter);
       else
	  printf("preconditioned Krylov subspace solver successfully completed with %ld iteration steps\n", iter);
    }

    for (j=0; j<A.nr; j++)
        printf("%12.4le",sol[j].r);
    printf("\n");
    for (j=0; j<A.nr; j++)
        printf("%12.4le",sol[j].i);
    printf("\n");fflush(stdout);


    
    //-------------------------------------------------------------------------
    // ------------------- clean up block triangular factors ------------------
    // ------------------------------------------------------------------------
    JanusDelete(&PREC);

    

    return (0);
} // end main

/*
    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Version:

	JANUS Block ILU R1.0.  

    Notice:

	Copyright (c) 2019 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/

	$Id: zhpdjanus.c 6263 2020-05-15 07:04:30Z bolle $
*/
