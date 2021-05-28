#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <janus.h>


int main(int argc, char **argv)
{
    // prescribe a non-real sparse matrix
    long int      ncol[8]={  1,    3,    5,    1,   1,   2,   4,   3},
                  rowind0[1]={0},
	          rowind1[3]=     {1,4,6},
		  rowind2[5]=           {0,1,2,5,7},
		  rowind3[1]=                 {3},
		  rowind4[1]=                      {1},
		  rowind5[2]=                           {0,5},
		  rowind6[4]=                                {0,3,6,7},
		  rowind7[3]=                                     {2,5,7},
                  *rowind[8]={rowind0,rowind1,rowind2,rowind3,rowind4,rowind5,rowind6,rowind7};
    doublecomplex val0[1]   ={{7.0,0.0}},
                  val1[3]   ={{-4.0,0.0},{-4.0,0.0},{1.0,0.0}},
		  val2[5]   ={{1.0,-1.0},{8.0,-1.0},{1.0,0.0},{7.0,-1.0},{-3.0,-1.0}},
		  val3[1]   ={{7.0,0.0}},
		  val4[1]   ={{2.0,-1.0}},
		  val5[2]   ={{2.0,-1.0},{3.0,-1.0}},
		  val6[4]   ={{7.0,-1.0},{9.0,-1.0},{11.0,-1.0},{2.0,-1.0}},
		  val7[3]   ={{5.0,-1.0},{8.0,-1.0},{5.0,-1.0}},
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
    A.nnz=20;
    A.ncol=ncol;
    A.rowind=rowind;
    A.val=(void *)val;
    /* tag properties of A */
    A.isreal=0;
    A.issymmetric=0;
    A.ishermitian=0;
    A.isdefinite=0;
    
    PrintMatrix(&A);



    //-------------------------------------------------------------------------
    //-------------- set up default user options for block ILU ---------------
    //-------------------------------------------------------------------------
    JanusDefaultOptions(&options);



    //-------------------------------------------------------------------------
    //---------------------------- compute block ILU -------------------------
    //-------------------------------------------------------------------------
    ierr=JanusFactor(&A, &PREC, options);
    if (ierr) {
       printf("computation of block incomplete CHOL failed\n");
       JanusDelete(&PREC);
       return (-1);
    }


    
    //-------------------------------------------------------------------------
    // ----- Krylov subspace solver (here: GMRES) preconditioned by BILU ------
    // ------------------------------------------------------------------------
    ierr=JanusSolver(&A,&PREC,rhs,sol,30,1e-6,5000,&iter);
    // why did the Krylov subspace solver stop?
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

	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://bilu.tu-bs.de/

	$Id: zgnljanus.c 6263 2020-05-15 07:04:30Z bolle $
*/
