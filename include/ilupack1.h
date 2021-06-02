/*! \file ilupack.h
   \brief main header for ILUPACK

   This header contains all definitions of functions as well as those of the
   constants
   $Id: ilupack1.h 826 2016-08-31 12:40:36Z bolle $
*/
#ifndef _ILU_PACK_H
#define _ILU_PACK_H
#define USE_LAPACK_DRIVER


#include <stdlib.h>

#include "long_integer.h"
#include "namesilupack.h"



#define LONG_INT integer
#define MEDIUM_INT int



#define PERM_NONE        0
#define PERM_MTMETIS     1
#define PERM_METISN      2
#define PERM_METISE      3
#define PERM_AMD         4
#define PERM_RCM         5

#define BLOCK_NONE           0
#define BLOCK_ILU1T          1
#define BLOCK_ILUPT          2
#define BLOCK_SUPERNODES     3
#define BLOCK_APP_SUPERNODES 4
/* #define MIN_BLOCK_SIZE_JANUS 4 */
#define MIN_BLOCK_SIZE_JANUS 1 /* */
#define MAX_BLOCK_SIZE_JANUS 1024
/* #define MAX_BLOCK_SIZE_JANUS 65536 */
/* #define GROW_BLOCK_SIZE_JANUS 1.2 */
#define GROW_BLOCK_SIZE_JANUS 1.0


/* sparse matrix data structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;             /* total number of columns */
  integer nnz;            /* total number of nonzeros */
  integer *ncol;          /* number of elements per column */
  integer **rowind;       /* row indices for each column */
  doubleprecision **val;  /* numerical values for each column */
} DSparseMatrix;
/* sparse lower triangular block matrix with one diagonal block and one lower 
   off-diagonal block per block column 
*/
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;		  /* total number of columns */
  integer nnz;            /* number of nonzeros */
  integer nblocks;	  /* total number of blocks */
  integer *nblockcol;     /* number columns per diagonal block */
  integer *nblockrow;     /* number of rows per sub-diagonal block */
  integer **colind;       /* column indices for each diagonal block */
  integer **rowind;       /* indices for the sub-diagonal block entries */
  doubleprecision **valD; /* numerical values of the square diagonal blocks */
  doubleprecision **valE; /* numerical values of the sub-diagonal blocks */
} DSparseBlockMatrix;


/* sparse matrix data structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;             /* total number of columns */
  integer nnz;            /* total number of nonzeros */
  integer *ncol;          /* number of elements per column */
  integer **rowind;       /* row indices for each column */
  doublecomplex **val;    /* numerical values for each column */
} ZSparseMatrix;
/* sparse lower triangular block matrix with one diagonal block and one lower
   off-diagonal block per block column
*/
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows    */
  integer nc;		  /* total number of columns */
  integer nnz;            /* number of nonzeros */
  integer nblocks;	  /* total number of blocks  */
  integer *nblockcol;     /* number columns per diagonal block */
  integer *nblockrow;     /* number of rows per sub-diagonal block */
  integer **colind;       /* column indices for each diagonal block */
  integer **rowind;       /* indices for the sub-diagonal block entries */
  doublecomplex **valD;   /* numerical values of the square diagonal blocks */
  doublecomplex **valE;   /* numerical values of the sub-diagonal blocks */
} ZSparseBlockMatrix;


/* sparse matrix data structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;             /* total number of columns */
  integer nnz;            /* total number of nonzeros */
  integer *ncol;          /* number of elements per column */
  integer **rowind;       /* row indices for each column */
  real    **val;          /* numerical values for each column */
} SSparseMatrix;
/* sparse lower triangular block matrix with one diagonal block and one lower
   off-diagonal block per block column
*/
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows   */
  integer nc;		  /* total number of columns */
  integer nnz;            /* number of nonzeros */
  integer nblocks;	  /* total number of blocks */
  integer *nblockcol;     /* number columns per diagonal block */
  integer *nblockrow;     /* number of rows per sub-diagonal block */
  integer **colind;       /* column indices for each diagonal block */
  integer **rowind;       /* indices for the sub-diagonal block entries */
  real **valD;            /* numerical values of the square diagonal blocks */
  real **valE;            /* numerical values of the sub-diagonal blocks */
} SSparseBlockMatrix;


/* sparse matrix data structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;             /* total number of columns */
  integer nnz;            /* total number of nonzeros */
  integer *ncol;          /* number of elements per column */
  integer **rowind;       /* row indices for each column */
  complex **val;          /* numerical values for each column */
} CSparseMatrix;
/* sparse lower triangular block matrix with one diagonal block and one lower
   off-diagonal block per block column
*/
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows    */
  integer nc;		  /* total number of columns */
  integer nnz;            /* number of nonzeros */
  integer nblocks;	  /* total number of blocks */
  integer *nblockcol;     /* number columns per diagonal block */
  integer *nblockrow;     /* number of rows per sub-diagonal block */
  integer **colind;       /* column indices for each diagonal block */
  integer **rowind;       /* indices for the sub-diagonal block entries */
  complex **valD;         /* numerical values of the square diagonal blocks */
  complex **valE;         /* numerical values of the sub-diagonal blocks */
} CSparseBlockMatrix;


/* sparse matrix data structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
  integer nr;             /* total number of rows */
  integer nc;             /* total number of columns */
  integer nnz;            /* total number of nonzeros */
  integer *ncol;          /* number of elements per column */
  integer **rowind;       /* row indices for each column */
  void **val;             /* numerical values for each column */
} SparseMatrix;

/* sparse lower triangular block matrix with one diagonal block and one lower
   off-diagonal block per block column
*/
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;             /* total number of rows    */
   integer nc;		  /* total number of columns */
   integer nnz;            /* number of nonzeros */
   integer nblocks;	  /* total number of blocks */
   integer *nblockcol;     /* number columns per diagonal block */
   integer *nblockrow;     /* number of rows per sub-diagonal block */
   integer **colind;       /* column indices for each diagonal block */
   integer **rowind;       /* indices for the sub-diagonal block entries */
   void **valD;         /* numerical values of the square diagonal blocks */
   void **valE;         /* numerical values of the sub-diagonal blocks */
} SparseBlockMatrix;


/* BILU options */
typedef struct {
  integer matching;       /* turn matching on/off */
  integer ordering;       /* which reorodering to take */
  double  droptol;	  /* drop tolerance  */
  integer cosine;         /* cosine-based blocking turned on/off */
  integer blocking_strategy; /* blocking strategy: ILU(1,droptol), supernodes,... */
  integer progressive_aggregation;/* progressively aggreate blocks during the factorization on/off */
  integer *blocksize;     /* optionally, pass a fixed block partitioning*/
  integer nblocks;        /* number of blocks */
  integer perturbation;   /* perturb singular diagonal blocks artificially */
  integer symmetric_structure; /* use symmetric permutations when the nonsymmetric matrix has almost symmetric structure */
  integer invert_blocks; /* invert diagonal blocks explicitly */
  integer level_of_fill; /* level-of-fill parameter for ILUPT */
} JanusOptions;

/* complete BILU preconditioning structure */
typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer n;
   integer invert_blocks; 
   doublecomplex logdet;
  
   SparseBlockMatrix *BL;
   SparseBlockMatrix *BiD;
   SparseBlockMatrix *BUT;
  
   void *SL;
   void *SR;
  
   integer *p;
   integer *invq;
   integer *pivots; 

} JanusPrec;




#define DGNLbilu  Dbilu
#define ZGNLbilu  Zbilu
#define DGNLbilua Dbilua
#define ZGNLbilua Zbilua

_CPP_PREFIX void JanusDefaultOptions(JanusOptions *options);
_CPP_PREFIX void JanusDelete(JanusPrec *);
_CPP_PREFIX void JanusInit(JanusPrec *);
_CPP_PREFIX integer JanusFactor(SparseMatrix *A, JanusPrec *PREC, JanusOptions options);
_CPP_PREFIX void JanusSol(JanusPrec *PREC, void *x, void *y, void *buff, integer m);
_CPP_PREFIX void JanusSolH(JanusPrec *PREC, void *x, void *y, void *buff, integer m);
_CPP_PREFIX void JanusSolT(JanusPrec *PREC, void *x, void *y, void *buff, integer m);
_CPP_PREFIX void JanusNnz(JanusPrec *PREC, integer *nnz, integer *mxblock, double *avblock, double *stddev);
_CPP_PREFIX integer JanusSolver(SparseMatrix *A, JanusPrec *PREC, void *rhs, void *sol,
				integer K_RESTART, double RES_TOL, integer MAX_IT, integer *iter);

_CPP_PREFIX void PrintMatrix(SparseMatrix *A);
_CPP_PREFIX void DPrintMatrix(DSparseMatrix *A);
_CPP_PREFIX void ZPrintMatrix(ZSparseMatrix *A);


_CPP_PREFIX integer DGNLbilu_driver(DSparseMatrix *A,
				    DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    JanusOptions options);
_CPP_PREFIX integer ZGNLbilu_driver(ZSparseMatrix *A,
				    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    JanusOptions options);
_CPP_PREFIX integer DSPDbilu_driver(DSparseMatrix *A,
				    DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    JanusOptions options);
_CPP_PREFIX integer DSYMbilu_driver(DSparseMatrix *A,
				    DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    integer *isdefinite,
				    JanusOptions options);
_CPP_PREFIX integer DSYMilu_driver(DSparseMatrix *A,
				   DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
				   double *SL, double *SR, integer *p, integer *invq,
				   doublecomplex *logdet,
				   integer *isdefinite,
				   JanusOptions options);
_CPP_PREFIX integer ZHPDbilu_driver(ZSparseMatrix *A,
				    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    JanusOptions options);
_CPP_PREFIX integer ZHERbilu_driver(ZSparseMatrix *A,
				    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    integer *isdefinite,
				    JanusOptions options);
_CPP_PREFIX integer ZHERilu_driver(ZSparseMatrix *A,
				   ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT,
				   double *SL, double *SR, integer *p, integer *invq,
				   doublecomplex *logdet,
				   integer *isdefinite,
				   JanusOptions options);
_CPP_PREFIX integer ZSYMbilu_driver(ZSparseMatrix *A,
				    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
				    double *SL, double *SR, integer *p, integer *invq,
				    doublecomplex *logdet,
				    integer *isdefinite,
				    JanusOptions options);
_CPP_PREFIX integer ZSYMilu_driver(ZSparseMatrix *A,
				   ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT,
				   double *SL, double *SR, integer *p, integer *invq,
				   doublecomplex *logdet,
				   integer *isdefinite,
				   JanusOptions options);
_CPP_PREFIX integer Dbilu(DSparseMatrix *A,
			  DSparseBlockMatrix *BL,
			  DSparseBlockMatrix *BiD,
			  DSparseBlockMatrix *BUT,
			  double *SL, double *SR, integer *p, integer *invq,
			  doublecomplex *logdet,
			  integer *blocksize, integer nblocks,
			  double droptol, integer ilu1t,
			  integer pa, integer perturbation,
			  integer symmetric_structure,
			  integer *pivots, integer level_of_fill);
_CPP_PREFIX integer Zbilu(ZSparseMatrix *A,
			  ZSparseBlockMatrix *BL,
			  ZSparseBlockMatrix *BiD,
			  ZSparseBlockMatrix *BUT,
			  double *SL, double *SR, integer *p, integer *invq,
			  doublecomplex *logdet,
			  integer *blocksize, integer nblocks,
			  double droptol, integer ilu1t,
			  integer pa, integer perturbation,
			  integer symmetric_structure,
			  integer *pivots, integer level_of_fill);
_CPP_PREFIX integer Dbilua(DSparseMatrix *A,
			   DSparseBlockMatrix *BL,
			   DSparseBlockMatrix *BiD,
			   DSparseBlockMatrix *BUT,
			   double *SL, double *SR, integer *p, integer *invq,
			   doublecomplex *logdet,
			   integer *blocksize, integer nblocks,
			   double droptol, integer ilu1t, integer pa,
			   integer perturbation,
			   integer symmetric_structure,
			   integer *pivots, integer level_of_fill);
_CPP_PREFIX integer Zbilua(ZSparseMatrix *A,
			   ZSparseBlockMatrix *BL,
			   ZSparseBlockMatrix *BiD,
			   ZSparseBlockMatrix *BUT,
			   double *SL, double *SR, integer *p, integer *invq,
			   doublecomplex *logdet,
			   integer *blocksize, integer nblocks,
			   double droptol, integer ilu1t,
			   integer pa, integer perturbation,
			   integer symmetric_structure,
			   integer *pivots, integer level_of_fill);
_CPP_PREFIX integer DSPDbilu(DSparseMatrix *A,
			     DSparseBlockMatrix *BL,
			     DSparseBlockMatrix *BiD,
			     DSparseBlockMatrix *BUT,
			     double *SL, double *SR, integer *p, integer *invq,
			     doublecomplex *logdet,
			     integer *blocksize, integer nblocks,
			     double droptol, integer ildl1t,
			     integer pa, integer perturbation,
			     integer invert_blocks, integer level_of_fill);
_CPP_PREFIX integer ZHPDbilu(ZSparseMatrix *A,
			     ZSparseBlockMatrix *BL,
			     ZSparseBlockMatrix *BiD,
			     ZSparseBlockMatrix *BUT,
			     double *SL, double *SR, integer *p, integer *invq,
			     doublecomplex *logdet,
			     integer *blocksize, integer nblocks,
			     double droptol, integer ildl1t,
			     integer pa, integer perturbation,
			     integer invert_blocks, integer level_of_fill);
_CPP_PREFIX integer DSYMbilu(DSparseMatrix *A,
			     DSparseBlockMatrix *BL,
			     DSparseBlockMatrix *BiD,
			     DSparseBlockMatrix *BUT,
			     double *SL, double *SR, integer *p, integer *invq,
			     doublecomplex *logdet,
			     integer *isdefinite,
			     integer *blocksize, integer nblocks,
			     double droptol, integer ildl1t,
			     integer pa, integer perturbation,
			     integer *pivots, integer level_of_fill);
_CPP_PREFIX integer ZHERbilu(ZSparseMatrix *A,
			     ZSparseBlockMatrix *BL,
			     ZSparseBlockMatrix *BiD,
			     ZSparseBlockMatrix *BUT,
			     double *SL, double *SR, integer *p, integer *invq,
			     doublecomplex *logdet,
			     integer *isdefinite,
			     integer *blocksize, integer nblocks,
			     double droptol, integer ildl1t,
			     integer pa, integer perturbation,
			     integer *pivots, integer level_of_fill);
_CPP_PREFIX integer ZSYMbilu(ZSparseMatrix *A,
			     ZSparseBlockMatrix *BL,
			     ZSparseBlockMatrix *BiD,
			     ZSparseBlockMatrix *BUT,
			     double *SL, double *SR, integer *p, integer *invq,
			     doublecomplex *logdet,
			     integer *isdefinite,
			     integer *blocksize, integer nblocks,
			     double droptol, integer ildl1t,
			     integer pa, integer perturbation,
			     integer *pivots, integer level_of_fill);
_CPP_PREFIX integer DSYMilu(DSparseMatrix *A,
			    DSparseBlockMatrix *BL,
			    DSparseBlockMatrix *BiD,
			    DSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq, 
			    doublecomplex *logdet,
			    integer *isdefinite,
			    double droptol);
_CPP_PREFIX integer ZHERilu(ZSparseMatrix *A,
			    ZSparseBlockMatrix *BL,
			    ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq, 
			    doublecomplex *logdet,
			    integer *isdefinite,
			    double droptol);
_CPP_PREFIX integer ZSYMilu(ZSparseMatrix *A,
			    ZSparseBlockMatrix *BL,
			    ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq, 
			    doublecomplex *logdet,
			    integer *isdefinite,
			    double droptol);
_CPP_PREFIX integer DSYMbilua(DSparseMatrix *A,
			      DSparseBlockMatrix *BL,
			      DSparseBlockMatrix *BiD,
			      DSparseBlockMatrix *BUT,
			      double *SL, double *SR, integer *p, integer *invq,
			      doublecomplex *logdet,
			      integer *isdefinite,
			      integer *blocksize, integer nblocks,
			      double droptol, integer ildl1t,
			      integer pa, integer perturbation,
			      integer *pivots, integer level_of_fill);
_CPP_PREFIX integer ZHERbilua(ZSparseMatrix *A,
			      ZSparseBlockMatrix *BL,
			      ZSparseBlockMatrix *BiD,
			      ZSparseBlockMatrix *BUT,
			      double *SL, double *SR, integer *p, integer *invq,
			      doublecomplex *logdet,
			      integer *isdefinite,
			      integer *blocksize, integer nblocks,
			      double droptol, integer ildl1t,
			      integer pa, integer perturbation,
			      integer *pivots, integer level_of_fill);
_CPP_PREFIX integer ZSYMbilua(ZSparseMatrix *A,
			      ZSparseBlockMatrix *BL,
			      ZSparseBlockMatrix *BiD,
			      ZSparseBlockMatrix *BUT,
			      double *SL, double *SR, integer *p, integer *invq,
			      doublecomplex *logdet,
			      integer *isdefinite,
			      integer *blocksize, integer nblocks,
			      double droptol, integer ildl1t,
			      integer pa, integer perturbation,
			      integer *pivots, integer level_of_fill);
_CPP_PREFIX void DSYMbfspai(DSparseMatrix *iA,
			    DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD,
			    DSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq,
			    double droptol);
_CPP_PREFIX void ZSYMbfspai(ZSparseMatrix *iA,
			    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq,
			    double droptol);
_CPP_PREFIX void ZHERbfspai(ZSparseMatrix *iA,
			    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq,
			    double droptol);

_CPP_PREFIX void DSYMselbinv(DSparseMatrix *iA,
			    DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD,
			    DSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq);
_CPP_PREFIX void ZSYMselbinv(ZSparseMatrix *iA,
			    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq);
_CPP_PREFIX void ZHERselbinv(ZSparseMatrix *iA,
			    ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD,
			    ZSparseBlockMatrix *BUT,
			    double *SL, double *SR, integer *p, integer *invq);


_CPP_PREFIX void Dbilusol(DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
			  double *rhssol, double *dbuff, integer n, integer m);
_CPP_PREFIX void Zbilusol(ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
			  doublecomplex *rhssol, doublecomplex *dbuff, integer n, integer m);
_CPP_PREFIX void Dbilusolt(DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
			   double *rhssol, double *dbuff, integer n, integer m);
_CPP_PREFIX void Zbilusolt(ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
			   doublecomplex *rhssol, doublecomplex *dbuff, integer n, integer m);
_CPP_PREFIX void Dbilusolh(DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
			   double *rhssol, double *dbuff, integer n, integer m);
_CPP_PREFIX void Zbilusolh(ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
			   doublecomplex *rhssol, doublecomplex *dbuff, integer n, integer m);
_CPP_PREFIX void DSYMbilusol(DSparseBlockMatrix *BL, DSparseBlockMatrix *BiD, DSparseBlockMatrix *BUT, integer *pivots,
			     double *rhssol, double *dbuff, integer n, integer m);
_CPP_PREFIX void ZSYMbilusol(ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
			     doublecomplex *rhssol, doublecomplex *dbuff, integer n, integer m);
_CPP_PREFIX void ZHERbilusol(ZSparseBlockMatrix *BL, ZSparseBlockMatrix *BiD, ZSparseBlockMatrix *BUT, integer *pivots,
			     doublecomplex *rhssol, doublecomplex *dbuff, integer n, integer m);
_CPP_PREFIX integer Dilu1t_blocks(DSparseMatrix *A,
				  double *SL, double *SR, integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double droptol);
_CPP_PREFIX integer Dilupt_blocks(DSparseMatrix *A,
				  double *SL, double *SR, integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double droptol, integer level_of_fill);
_CPP_PREFIX integer Zilupt_blocks(ZSparseMatrix *A,
				  double *SL, double *SR, integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double droptol, integer level_of_fill);
_CPP_PREFIX integer Dildlpt_blocks(DSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol, integer level_of_fill);
_CPP_PREFIX integer DSYMbildlpt_blocks(DSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol, integer level_of_fill);
_CPP_PREFIX integer Zildlpt_blocks(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol, integer level_of_fill);
_CPP_PREFIX integer ZSYMbildlpt_blocks(ZSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol, integer level_of_fill);
_CPP_PREFIX integer ZHERbildlpt_blocks(ZSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol, integer level_of_fill);
_CPP_PREFIX integer Dsupernodes(DSparseMatrix *A,
				double *SL, double *SR, integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double droptol);
_CPP_PREFIX integer DSPDappsupernodes(DSparseMatrix *A,
				      double *SL, double *SR, integer *p, integer *invq,
				      integer *blocksize, integer *nblocks,
				      double droptol);
_CPP_PREFIX integer DSYMsupernodes(DSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Zilu1t_blocks(ZSparseMatrix *A,
				  double *SL, double *SR, integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double droptol);
_CPP_PREFIX integer Zsupernodes(ZSparseMatrix *A,
				double *SL, double *SR, integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double droptol);
_CPP_PREFIX integer ZHPDappsupernodes(ZSparseMatrix *A,
				      double *SL, double *SR, integer *p, integer *invq,
				      integer *blocksize, integer *nblocks,
				      double droptol);
_CPP_PREFIX integer ZSYMsupernodes(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer ZHERsupernodes(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Dilu_blocks(DSparseMatrix *A,
				double *SL, double *SR, integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double droptol);
_CPP_PREFIX integer Zilu_blocks(ZSparseMatrix *A,
				double *SL, double *SR, integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double droptol);
_CPP_PREFIX integer Dbilu1t_blocks(DSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Zbilu1t_blocks(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Dildl1t_blocks(DSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Zildl1t_blocks(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer DSYMbildl1t_blocks(DSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol);
_CPP_PREFIX integer ZSYMbildl1t_blocks(ZSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol);
_CPP_PREFIX integer ZHERbildl1t_blocks(ZSparseMatrix *A,
				       double *SL, double *SR, integer *p, integer *invq,
				       integer *blocksize, integer *nblocks,
				       double droptol);
_CPP_PREFIX integer DSYMildl_blocks(DSparseMatrix *A,
				    double *SL, double *SR, integer *p, integer *invq,
				    integer *blocksize, integer *nblocks,
				    double droptol);
_CPP_PREFIX integer ZSYMildl_blocks(ZSparseMatrix *A,
				    double *SL, double *SR, integer *p, integer *invq,
				    integer *blocksize, integer *nblocks,
				    double droptol);
_CPP_PREFIX integer ZHERildl_blocks(ZSparseMatrix *A,
				    double *SL, double *SR, integer *p, integer *invq,
				    integer *blocksize, integer *nblocks,
				    double droptol);
_CPP_PREFIX integer DSYMbildl1t_blocksa(DSparseMatrix *A,
					double *SL, double *SR, integer *p, integer *invq,
					integer *blocksize, integer *nblocks,
					double droptol);
_CPP_PREFIX integer ZSYMbildl1t_blocksa(ZSparseMatrix *A,
					double *SL, double *SR, integer *p, integer *invq,
					integer *blocksize, integer *nblocks,
					double droptol);
_CPP_PREFIX integer ZHERbildl1t_blocksa(ZSparseMatrix *A,
					double *SL, double *SR, integer *p, integer *invq,
					integer *blocksize, integer *nblocks,
					double droptol);
_CPP_PREFIX integer Dilu1t_blocksa(DSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX integer Zilu1t_blocksa(ZSparseMatrix *A,
				   double *SL, double *SR, integer *p, integer *invq,
				   integer *blocksize, integer *nblocks,
				   double droptol);
_CPP_PREFIX void Dcosine_blocks(DSparseMatrix *A,
				integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double tau);
_CPP_PREFIX void Zcosine_blocks(ZSparseMatrix *A,
				integer *p, integer *invq,
				integer *blocksize, integer *nblocks,
				double tau);
_CPP_PREFIX void Dmcosine_blocks(DSparseMatrix *A,
				 integer *p, integer *invq,
				 integer *blocksize, integer *nblocks,
				 double tau);
_CPP_PREFIX void Zmcosine_blocks(ZSparseMatrix *A,
				 integer *p, integer *invq,
				 integer *blocksize, integer *nblocks,
				 double tau);
_CPP_PREFIX void Dmcosine_sblocks(DSparseMatrix *A,
				  integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double tau);
_CPP_PREFIX void Zmcosine_sblocks(ZSparseMatrix *A,
				  integer *p, integer *invq,
				  integer *blocksize, integer *nblocks,
				  double tau);


_CPP_PREFIX void MatVec(SparseMatrix *A, void *x, void *y);
_CPP_PREFIX void DMatVec(DSparseMatrix *A, double *x, double *y);
_CPP_PREFIX void ZMatVec(ZSparseMatrix *A, doublecomplex *x, doublecomplex *y);
_CPP_PREFIX void DMatTVec(DSparseMatrix *A, double *x, double *y);
_CPP_PREFIX void ZMatTVec(ZSparseMatrix *A, doublecomplex *x, doublecomplex *y);
_CPP_PREFIX void DMatHVec(DSparseMatrix *A, double *x, double *y);
_CPP_PREFIX void ZMatHVec(ZSparseMatrix *A, doublecomplex *x, doublecomplex *y);
_CPP_PREFIX void DSYMMatVec(DSparseMatrix *A, double *x, double *y);
_CPP_PREFIX void ZSYMMatVec(ZSparseMatrix *A, doublecomplex *x, doublecomplex *y);
_CPP_PREFIX void ZHERMatVec(ZSparseMatrix *A, doublecomplex *x, doublecomplex *y);

_CPP_PREFIX void SparseBlockDelete(SparseBlockMatrix *B);
_CPP_PREFIX void DSparseBlockDelete(DSparseBlockMatrix *B);
_CPP_PREFIX void ZSparseBlockDelete(ZSparseBlockMatrix *B);

_CPP_PREFIX void SparseBlockInit(SparseBlockMatrix *B);
_CPP_PREFIX void DSparseBlockInit(DSparseBlockMatrix *B);
_CPP_PREFIX void ZSparseBlockInit(ZSparseBlockMatrix *B);

_CPP_PREFIX void SparseDelete(SparseMatrix *B);
_CPP_PREFIX void DSparseDelete(DSparseMatrix *B);
_CPP_PREFIX void ZSparseDelete(ZSparseMatrix *B);

_CPP_PREFIX void SparseInit(SparseMatrix *B);
_CPP_PREFIX void DSparseInit(DSparseMatrix *B);
_CPP_PREFIX void ZSparseInit(ZSparseMatrix *B);







/*! switch to indicate inverse-based dropping. It is used in AMGINIT
    AMGGETPARAMS and AMGSETPARAMS. The parameter "flag" is bitwise
    modified by flag|=DROP_INVERSE to set and flag&=~DROP_INVERSE to
    turn off inverse-based dropping.
    In AMGINIT, DROP_INVERSE is set by default
 */
#define DROP_INVERSE                     1

/*! switch for not shifting away zero pivots. This switch is used in ILUC
    which does not have pivoting prevent small diagonal entries.
    The parameter "param" is bitwise modified by param|=NO_SHIFT to 
    suppress shifts and param&=~NO_SHIFT to allow shifts.
 */
#define NO_SHIFT                         2

/*! switch for using Tismenetsky update. 
 */
#define TISMENETSKY_SC                   4
/* switch for repeated ILU */
#define REPEAT_FACT                      8
/* switch for enhanced estimate for the norm of the inverses */
#define IMPROVED_ESTIMATE               16
/* switch for using diagonal compensation */
#define DIAGONAL_COMPENSATION           32
/* switch for reducing the partial factorization to the non-coarse part */
#define COARSE_REDUCE                   64

/* switch for using a different pivoting strategy, if the regular reordering
   fails and before we switch to ILUTP
*/
#define FINAL_PIVOTING                 128
/* enforce the positve definite property */
#define ENSURE_SPD                     256

/* switch for the most simple Schur complement update */
#define SIMPLE_SC                      512


#define PREPROCESS_INITIAL_SYSTEM     1024
#define PREPROCESS_SUBSYSTEMS         2048
#define MULTI_PILUC                   4096


#define RE_FACTOR                     8192
#define AGGRESSIVE_DROPPING          16384
#define DISCARD_MATRIX               32768
#define SYMMETRIC_STRUCTURE          65536

#define STATIC_TESTVECTOR           131072
#define DYNAMIC_TESTVECTOR          262144
#define ADAPT_CONDEST               524288

#define SADDLE_POINT               1048576
#define BLOCK_STRUCTURE            2097152
#define DECOUPLE_CONSTRAINTS       4194304
#define DECOUPLE_CONSTRAINTSHH     8388608
#define NO_INITIAL_PREPROCESSING  16777216

/*                                 
                                  
                                  33554432
                                  67108864
                                 134217728 
                                 268435456 
                                 536870912
                                1073741824
                                2147483648
*/




#define _D_REAL_MAX_        1.7e+308
#define _S_REAL_MAX_        1.7e+38


/* ***************************************************** */
/* ******      Definitions for preconditioners     ***** */
typedef struct {
  integer *tree;
  integer *rang;
  integer *size;
  integer *chld;
  integer *nmch;
  integer *brth;
  integer *wght;
  integer *mark;
  integer *ownr;
  integer *hght;
  integer dimL;
  integer dimT;
  integer *perm;
  integer dimperm;
} OMPtab;

typedef struct {
  integer *tree;
  integer *rang;
  integer *size;
  integer *chld;
  integer *nmch;
  integer *brth;
  integer *wght;
  integer *mark;
  integer *ownr;
  integer *hght;
  integer dimL;
  integer dimT;
  integer *perm;
  integer dimperm;
} SOMPtab;

typedef struct {
  integer *tree;
  integer *rang;
  integer *size;
  integer *chld;
  integer *nmch;
  integer *brth;
  integer *wght;
  integer *mark;
  integer *ownr;
  integer *hght;
  integer dimL;
  integer dimT;
  integer *perm;
  integer dimperm;
} DOMPtab;

typedef struct {
  integer *tree;
  integer *rang;
  integer *size;
  integer *chld;
  integer *nmch;
  integer *brth;
  integer *wght;
  integer *mark;
  integer *ownr;
  integer *hght;
  integer dimL;
  integer dimT;
  integer *perm;
  integer dimperm;
} COMPtab;

typedef struct {
  integer *tree;
  integer *rang;
  integer *size;
  integer *chld;
  integer *nmch;
  integer *brth;
  integer *wght;
  integer *mark;
  integer *ownr;
  integer *hght;
  integer dimL;
  integer dimT;
  integer *perm;
  integer dimperm;
} ZOMPtab;

typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   void *a;
} SPARSEmat;

typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   doubleprecision *a;
} Dmat;

typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   real *a;
} Smat;

typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   doublecomplex *a;
} Zmat;

typedef struct {
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   complex *a;
} Cmat;



#define ILUPACK_NIPAR   50
#define ILUPACK_NFPAR   50




typedef struct  AMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer ispartial;
  integer nlev;                  
  integer n;                  
  integer nB;
  SPARSEmat A; 
  SPARSEmat LU;
  integer *LUperm;
  SPARSEmat E;
  SPARSEmat F;
  integer *p;
  integer *invq;
  void *rowscal;
  void *colscal;
  void *absdiag;
  struct AMGLM *prev;
  struct AMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  doubleprecision **blockdiag;
  integer **blockdiagind;
  integer *blockdiagsize;
  integer maxblockdiagsize;
  integer nblockdiag;
  integer **indB;
  doubleprecision  **LB;
  integer **indB2;
  doubleprecision   **LB2;
  integer *cntB;
  integer cnt;
  integer type_decouple;
  
  struct AMGLM *ompparts;
  struct AMGLM *omproot;
  integer nompparts;
  OMPtab omptab;
  integer *p_local;
  /*
     size_t *pardiso_pt[64];
     int     pardiso_iparm[64];
  */
} AMGlevelmat; 

typedef struct  DAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer ispartial;
  integer nlev;                  
  integer n;                  
  integer nB;
  Dmat A; 
  Dmat LU;
  integer *LUperm;
  Dmat E;
  Dmat F;
  integer *p;
  integer *invq;
  doubleprecision *rowscal;
  doubleprecision *colscal;
  doubleprecision *absdiag;
  struct DAMGLM *prev;
  struct DAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  doubleprecision **blockdiag;
  integer **blockdiagind;
  integer *blockdiagsize;
  integer maxblockdiagsize;
  integer nblockdiag;
  integer **indB;
  doubleprecision  **LB;
  integer **indB2;
  doubleprecision   **LB2;
  integer *cntB;
  integer cnt;
  integer type_decouple;

  struct DAMGLM *ompparts;
  struct DAMGLM *omproot;
  integer nompparts;
  DOMPtab omptab;
  integer *p_local;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} DAMGlevelmat; 

typedef struct  SAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer ispartial;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Smat A; 
  Smat LU;
  integer *LUperm;
  Smat E;
  Smat F;
  integer *p;
  integer *invq;
  real *rowscal;
  real *colscal;
  real *absdiag;
  struct SAMGLM *prev;
  struct SAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  real errorL;
  real errorU;
  real errorS;
  real **blockdiag;
  integer **blockdiagind;
  integer *blockdiagsize;
  integer maxblockdiagsize;
  integer nblockdiag;
  integer **indB;
  real  **LB;
  integer **indB2;
  real   **LB2;
  integer *cntB;
  integer cnt;
  integer type_decouple;

  struct SAMGLM *ompparts;
  struct SAMGLM *omproot;
  integer nompparts;
  SOMPtab omptab;
  integer *p_local;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} SAMGlevelmat; 

typedef struct  ZAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer ispartial;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Zmat A; 
  Zmat LU;
  integer *LUperm;
  Zmat E;
  Zmat F;
  integer *p;
  integer *invq;
  doublecomplex *rowscal;
  doublecomplex *colscal;
  doublecomplex *absdiag;
  struct ZAMGLM *prev;
  struct ZAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  doublecomplex **blockdiag;
  integer **blockdiagind;
  integer *blockdiagsize;
  integer maxblockdiagsize;
  integer nblockdiag;
  integer **indB;
  doublecomplex  **LB;
  integer **indB2;
  doublecomplex   **LB2;
  integer *cntB;
  integer cnt;
  integer type_decouple;

  struct ZAMGLM *ompparts;
  struct ZAMGLM *omproot;
  integer nompparts;
  ZOMPtab omptab;
  integer *p_local;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} ZAMGlevelmat; 

typedef struct CAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer ispartial;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Cmat A; 
  Cmat LU;
  integer *LUperm;
  Cmat E;
  Cmat F;
  integer *p;
  integer *invq;
  complex *rowscal;
  complex *colscal;
  complex *absdiag;
  struct CAMGLM *prev;
  struct CAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  real errorL;
  real errorU;
  real errorS;
  complex **blockdiag;
  integer **blockdiagind;
  integer *blockdiagsize;
  integer maxblockdiagsize;
  integer nblockdiag;
  integer **indB;
  complex **LB;
  integer **indB2;
  complex **LB2;
  integer *cntB;
  integer cnt;
  integer type_decouple;

  struct CAMGLM *ompparts;
  struct CAMGLM *omproot;
  integer nompparts;
  COMPtab omptab;
  integer *p_local;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} CAMGlevelmat; 


typedef struct {
   integer       ipar[ILUPACK_NIPAR];
   real          fpar[ILUPACK_NFPAR];
   integer       type;
   integer       *ibuff;
   integer       *iaux;
   real          *dbuff;
   real          *daux;
   integer       *ju;
   integer       *jlu;
   real          *alu;
   real          *testvector;
   size_t        nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer       rcomflag, returnlabel;
   real          *tv;
   integer       *ind;
   integer       nblocks;
   integer       *blocksize;
   integer       nindicator;
   integer       *indicator;
   Smat          A;
   integer       istack[30], *pistack[20];
   real          rstack[30], *prstack[10];
   real          fstack[30], *pfstack[10];
   size_t        ststack[5], *pststack[5];
   Smat          mstack[1];
   SAMGlevelmat  *amglmstack[1];
   integer       (*intfctstack[3])();
   integer       matching;
   char          *ordering;
   real          droptol;
   real          droptolS;
   real          droptolc;
   real          condest;
   real          restol;
   integer       maxit;
   real          elbow;
   integer       lfil;
   integer       lfilS;
   char          *typetv;
   char          *amg;
   integer       npresmoothing;
   integer       npostsmoothing;
   integer       ncoarse;
   char          *presmoother;
   char          *postsmoother;
   char          *FCpart;
   char          *typecoarse;
   integer       nrestart;
   integer       flags;
   char          *solver;
   real          damping;
   real          contraction;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   real          shift0;
   real          shiftmax;
   integer       nshifts;
   real          *shifts;
   Smat          *shiftmatrix;
   integer          niter;
   integer          nthreads;
   double           loadbalancefactor;
   integer          *indpartial;
   integer          *indexpartial;
   real             *valpartial;
   real             *valuepartial;
   real             *vpartial;
   real             *tvpartial;
} SILUPACKparam;

typedef struct {
   integer          ipar[ILUPACK_NIPAR];
   doubleprecision  fpar[ILUPACK_NFPAR];
   integer          type;
   integer          *ibuff;
   integer          *iaux;
   doubleprecision  *dbuff;
   doubleprecision  *daux;
   integer          *ju;
   integer          *jlu;
   doubleprecision  *alu;
   doubleprecision  *testvector;
   size_t           nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   doubleprecision  *tv;
   integer          *ind;
   integer       nblocks;
   integer       *blocksize;
   integer          nindicator;
   integer          *indicator;
   Dmat             A;
   integer          istack[30], *pistack[20];
   doubleprecision  rstack[30], *prstack[10];
   doubleprecision  fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Dmat             mstack[1];
   DAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer          matching;
   char             *ordering;
   doubleprecision  droptol;
   doubleprecision  droptolS;
   doubleprecision  droptolc;
   doubleprecision  condest;
   doubleprecision  restol;
   integer          maxit;
   doubleprecision  elbow;
   integer          lfil;
   integer          lfilS;
   char             *typetv;
   char             *amg;
   integer          npresmoothing;
   integer          npostsmoothing;
   integer          ncoarse;
   char             *presmoother;
   char             *postsmoother;
   char             *FCpart;
   char             *typecoarse;
   integer          nrestart;
   integer          flags;
   char             *solver;
   doubleprecision  damping;
   doubleprecision  contraction;
   integer          mixedprecision;
   integer          (*perm0)();
   integer          (*perm)();
   integer          (*permf)();
   doubleprecision  shift0;
   doubleprecision  shiftmax;
   integer          nshifts;
   doubleprecision  *shifts;
   Dmat             *shiftmatrix;
   integer          niter;
   integer          nthreads;
   double           loadbalancefactor;
   integer          *indpartial;
   integer          *indexpartial;
   doubleprecision  *valpartial;
   doubleprecision  *valuepartial;
   doubleprecision  *vpartial;
   doubleprecision  *tvpartial;
} DILUPACKparam;

typedef struct {
   integer     ipar[ILUPACK_NIPAR];
   real    fpar[ILUPACK_NFPAR];
   integer     type;
   integer     *ibuff;
   integer     *iaux;
   complex *dbuff;
   complex *daux;
   integer     *ju;
   integer     *jlu;
   complex *alu;
   complex *testvector;
   size_t  nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   complex          *tv;
   integer          *ind;
   integer       nblocks;
   integer       *blocksize;
   integer              nindicator;
   integer          *indicator;
   Cmat             A;
   integer          istack[30], *pistack[20];
   real             rstack[30], *prstack[10];
   complex          fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Cmat             mstack[1];
   CAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer      matching;
   char         *ordering;
   real         droptol;
   real         droptolS;
   real         droptolc;
   real         condest;
   real         restol;
   integer      maxit;
   real         elbow;
   integer      lfil;
   integer      lfilS;
   char         *typetv;
   char         *amg;
   integer      npresmoothing;
   integer      npostsmoothing;
   integer      ncoarse;
   char         *presmoother;
   char         *postsmoother;
   char         *FCpart;
   char         *typecoarse;
   integer       nrestart;
   integer       flags;
   char         *solver;
   real          damping;
   real          contraction;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   complex       shift0;
   complex       shiftmax;
   integer       nshifts;
   complex       *shifts;
   Cmat          *shiftmatrix;
   integer          niter;
   integer          nthreads;
   double           loadbalancefactor;
   integer          *indpartial;
   integer          *indexpartial;
   complex          *valpartial;
   complex          *valuepartial;
   complex          *vpartial;
   complex          *tvpartial;
} CILUPACKparam;

typedef struct {
   integer              ipar[ILUPACK_NIPAR];
   doubleprecision  fpar[ILUPACK_NFPAR];
   integer              type;
   integer              *ibuff;
   integer              *iaux;
   doublecomplex    *dbuff;
   doublecomplex    *daux;
   integer              *ju;
   integer              *jlu;
   doublecomplex    *alu;
   doublecomplex    *testvector;
   size_t           nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   doublecomplex    *tv;
   integer          *ind;
   integer       nblocks;
   integer       *blocksize;
   integer              nindicator;
   integer          *indicator;
   Zmat             A;
   integer          istack[30], *pistack[20];
   doubleprecision  rstack[30], *prstack[10];
   doublecomplex    fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Zmat             mstack[1];
   ZAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer      matching;
   char         *ordering;
   doubleprecision       droptol;
   doubleprecision       droptolS;
   doubleprecision       droptolc;
   doubleprecision       condest;
   doubleprecision       restol;
   integer      maxit;
   doubleprecision       elbow;
   integer      lfil;
   integer      lfilS;
   char         *typetv;
   char         *amg;
   integer      npresmoothing;
   integer      npostsmoothing;
   integer      ncoarse;
   char         *presmoother;
   char         *postsmoother;
   char         *FCpart;
   char         *typecoarse;
   integer       nrestart;
   integer       flags;
   char         *solver;
   doubleprecision damping;
   doubleprecision contraction;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   doublecomplex shift0;
   doublecomplex shiftmax;
   integer       nshifts;
   doublecomplex *shifts;
   Zmat          *shiftmatrix;
   integer          niter;
   integer          nthreads;
   double           loadbalancefactor;
   integer          *indpartial;
   integer          *indexpartial;
   doublecomplex    *valpartial;
   doublecomplex    *valuepartial;
   doublecomplex    *vpartial;
   doublecomplex    *tvpartial;
} ZILUPACKparam;


typedef struct {
   integer              ipar[ILUPACK_NIPAR];
   doubleprecision      fpar[ILUPACK_NFPAR];
   integer              type;
   integer              *ibuff;
   integer              *iaux;
   void                 *dbuff;
   void                 *daux;
   integer              *ju;
   integer              *jlu;
   void                 *alu;
   void                 *testvector;
   size_t               nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer              rcomflag, returnlabel;
   void                 *tv;
   integer              *ind;
   integer       nblocks;
   integer       *blocksize;
   integer              nindicator;
   integer              *indicator;
   SPARSEmat            A;
   integer              istack[30], *pistack[20];
   doubleprecision      rstack[30], *prstack[10];
   doublecomplex        fstack[30], *pfstack[10];
   size_t               ststack[5], *pststack[5];
   SPARSEmat            mstack[1];
   AMGlevelmat          *amglmstack[1];
   integer              (*intfctstack[3])();
   integer              matching;
   char                 *ordering;
   doubleprecision      droptol;
   doubleprecision      droptolS;
   doubleprecision      droptolc;
   doubleprecision      condest;
   doubleprecision      restol;
   integer              maxit;
   doubleprecision      elbow;
   integer              lfil;
   integer              lfilS;
   char                 *typetv;
   char                 *amg;
   integer              npresmoothing;
   integer              npostsmoothing;
   integer              ncoarse;
   char                 *presmoother;
   char                 *postsmoother;
   char                 *FCpart;
   char                 *typecoarse;
   integer              nrestart;
   integer              flags;
   char                 *solver;
   doubleprecision      damping;
   doubleprecision      contraction;
   integer              mixedprecision;
   integer              (*perm0)();
   integer              (*perm)();
   integer              (*permf)();
   integer              niter;
   integer          nthreads;
   double           loadbalancefactor;
   integer          *indpartial;
   integer          *indexpartial;
   void             *valpartial;
   void             *valuepartial;
   void             *vpartial;
   void             *tvpartial;
} ILUPACKparam;



_CPP_PREFIX void DGNLAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void SGNLAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
			     SILUPACKparam *, real *, real *);
_CPP_PREFIX void ZGNLAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CGNLAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
			     CILUPACKparam *, complex *, complex *);

_CPP_PREFIX void DGNLAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void SGNLAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
			     SILUPACKparam *, real *, real *);
_CPP_PREFIX void ZGNLAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CGNLAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
			     CILUPACKparam *, complex *, complex *);


_CPP_PREFIX void DSPDAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void SSPDAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     SILUPACKparam *, real *, real *);
_CPP_PREFIX void ZHPDAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CHPDAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *);

_CPP_PREFIX void DSPDAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void SSPDAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     SILUPACKparam *, real *, real *);
_CPP_PREFIX void ZHPDAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CHPDAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *);


_CPP_PREFIX void DSYMAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void SSYMAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     SILUPACKparam *, real *, real *);
_CPP_PREFIX void ZHERAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CHERAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void ZSYMAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void CSYMAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *);

_CPP_PREFIX void DSYMAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     DILUPACKparam *, doubleprecision *, doubleprecision *, 
			     integer *);
_CPP_PREFIX void SSYMAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     SILUPACKparam *, real *, real *, 
			     integer *);
_CPP_PREFIX void ZHERAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *, 
			     integer *);
_CPP_PREFIX void CHERAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *,
			     integer *);
_CPP_PREFIX void ZSYMAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     ZILUPACKparam *, doublecomplex *, doublecomplex *, 
			     integer *);
_CPP_PREFIX void CSYMAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
                             integer, integer, 
			     CILUPACKparam *, complex *, complex *, 
			     integer *);





_CPP_PREFIX void   AMGsol(AMGlevelmat *,ILUPACKparam *, void *, void *);
_CPP_PREFIX void   AMGtsol(AMGlevelmat *,ILUPACKparam *, void *, void *, void *);
_CPP_PREFIX void   AMGhsol(AMGlevelmat *,ILUPACKparam *, void *, void *, void *);


_CPP_PREFIX void    DGNLAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void    DGNLAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void  DGNLAMGdlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void  DGNLAMGdusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void   DGNLAMGlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void   DGNLAMGusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLAMGtdlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLAMGtdusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void  DGNLAMGtlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void  DGNLAMGtusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void   DGNLAMGtsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void   DGNLAMGtsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void    DSPDAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void    DSPDAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
void    DSYMAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void    DSYMAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void    DSYMAMGbsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLAMGextract(Dmat *,Dmat *, Dmat, integer *,integer *,  integer);
_CPP_PREFIX void DSYMAMGextract(Dmat *, Dmat, integer *,integer *,  integer);
_CPP_PREFIX void DSSMAMGextract(Dmat *, Dmat, integer *,integer *,  integer);

_CPP_PREFIX void    SGNLAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void    SGNLAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
_CPP_PREFIX void  SGNLAMGdlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void  SGNLAMGdusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void   SGNLAMGusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void   SGNLAMGlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void SGNLAMGtdlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void SGNLAMGtdusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void  SGNLAMGtusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void  SGNLAMGtlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void   SGNLAMGtsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void   SGNLAMGtsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
_CPP_PREFIX void    SSPDAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void    SSPDAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
_CPP_PREFIX void    SSYMAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void    SSYMAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
_CPP_PREFIX void    SSYMAMGbsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
_CPP_PREFIX void SGNLAMGextract(Smat *,Smat *, Smat, integer *,integer *,  integer);
_CPP_PREFIX void SSYMAMGextract(Smat *, Smat, integer *,integer *,  integer);
_CPP_PREFIX void SSSMAMGextract(Smat *, Smat, integer *,integer *,  integer);

_CPP_PREFIX void    ZGNLAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZGNLAMGsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void  ZGNLAMGdlsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void  ZGNLAMGdusol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void   ZGNLAMGlsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void   ZGNLAMGusol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZGNLAMGtdlsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZGNLAMGtdusol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void  ZGNLAMGtlsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void  ZGNLAMGtusol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void   ZGNLAMGtsol_internal(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void   ZGNLAMGtsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZHPDAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZHPDAMGsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZHERAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZHERAMGsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZSYMAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZSYMAMGsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZHERAMGbsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void    ZSYMAMGbsol(ZAMGlevelmat *,ZILUPACKparam *, doublecomplex *, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZGNLAMGextract(Zmat *,Zmat *, Zmat, integer *,integer *,  integer);
_CPP_PREFIX void ZHERAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
_CPP_PREFIX void ZSHRAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
_CPP_PREFIX void ZSYMAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
_CPP_PREFIX void ZSSMAMGextract(Zmat *, Zmat, integer *,integer *,  integer);

_CPP_PREFIX void    CGNLAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void    CGNLAMGsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void  CGNLAMGdlsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void  CGNLAMGdusol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void   CGNLAMGlsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void   CGNLAMGusol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void CGNLAMGtdlsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void CGNLAMGtdusol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void  CGNLAMGtusol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void  CGNLAMGtlsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void   CGNLAMGtsol_internal(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void   CGNLAMGtsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void    CHPDAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void    CHPDAMGsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void    CHERAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void    CHERAMGsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void    CSYMAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void    CSYMAMGsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *);
_CPP_PREFIX void    CHERAMGbsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void    CSYMAMGbsol(CAMGlevelmat *,CILUPACKparam *, complex *, complex *, complex *);
_CPP_PREFIX void CGNLAMGextract(Cmat *,Cmat *, Cmat, integer *,integer *,  integer);
_CPP_PREFIX void CHERAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
_CPP_PREFIX void CSHRAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
_CPP_PREFIX void CSYMAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
_CPP_PREFIX void CSSMAMGextract(Cmat *, Cmat, integer *,integer *,  integer);


_CPP_PREFIX void DGNLlupq(integer *, doubleprecision *, integer *, integer *, integer *);
_CPP_PREFIX void DGNLlupqsol (integer *, doubleprecision *, integer *, integer *, 
			      doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqtsol(integer *, doubleprecision *, integer *, integer *, 
			      doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqlsol (integer *, doubleprecision *, integer *, integer *, 
			       doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqtlsol(integer *, doubleprecision *, integer *, integer *, 
			       doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqusol (integer *, doubleprecision *, integer *, integer *, 
			       doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqtusol(integer *, doubleprecision *, integer *, integer *, 
			       doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqdlsol (integer *, doubleprecision *, integer *, integer *, 
				doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqtdlsol(integer *, doubleprecision *, integer *, integer *, 
				doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqdusol (integer *, doubleprecision *, integer *, integer *, 
				doubleprecision *, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLlupqtdusol(integer *, doubleprecision *, integer *, integer *, 
				doubleprecision *, doubleprecision *, doubleprecision *);

_CPP_PREFIX void SGNLlupq(integer *, real *, integer *, integer *, integer *);
_CPP_PREFIX void SGNLlupqsol (integer *, real *, integer *, integer *, real *, real *, 
			      real *);
_CPP_PREFIX void SGNLlupqtsol(integer *, real *, integer *, integer *, real *, real *, 
			      real *);
_CPP_PREFIX void SGNLlupqlsol (integer *, real *, integer *, integer *, real *, real *, 
			       real *);
_CPP_PREFIX void SGNLlupqtlsol(integer *, real *, integer *, integer *, real *, real *, 
			       real *);
_CPP_PREFIX void SGNLlupqusol (integer *, real *, integer *, integer *, real *, real *, 
			       real *);
_CPP_PREFIX void SGNLlupqtusol(integer *, real *, integer *, integer *, real *, real *, 
			       real *);
_CPP_PREFIX void SGNLlupqdlsol (integer *, real *, integer *, integer *, real *, real *, 
				real *);
_CPP_PREFIX void SGNLlupqtdlsol(integer *, real *, integer *, integer *, real *, real *, 
				real *);
_CPP_PREFIX void SGNLlupqdusol (integer *, real *, integer *, integer *, real *, real *, 
				real *);
_CPP_PREFIX void SGNLlupqtdusol(integer *, real *, integer *, integer *, real *, real *, 
				real *);

_CPP_PREFIX void ZGNLlupq(integer *, doublecomplex *, integer *, integer *, integer *);
_CPP_PREFIX void ZGNLlupqsol (integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			      doublecomplex *);
_CPP_PREFIX void ZGNLlupqtsol(integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			      doublecomplex *);
_CPP_PREFIX void ZGNLlupqlsol (integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			       doublecomplex *);
_CPP_PREFIX void ZGNLlupqtlsol(integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			       doublecomplex *);
_CPP_PREFIX void ZGNLlupqusol (integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			       doublecomplex *);
_CPP_PREFIX void ZGNLlupqtusol(integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
			       doublecomplex *);
_CPP_PREFIX void ZGNLlupqdlsol (integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
				doublecomplex *);
_CPP_PREFIX void ZGNLlupqtdlsol(integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
				doublecomplex *);
_CPP_PREFIX void ZGNLlupqdusol (integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
				doublecomplex *);
_CPP_PREFIX void ZGNLlupqtdusol(integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, 
				doublecomplex *);

_CPP_PREFIX void CGNLlupq(integer *, complex *, integer *, integer *, integer *);
_CPP_PREFIX void CGNLlupqsol (integer *, complex *, integer *, integer *, complex *, complex *, 
			      complex *);
_CPP_PREFIX void CGNLlupqtsol(integer *, complex *, integer *, integer *, complex *, complex *, 
			      complex *);
_CPP_PREFIX void CGNLlupqlsol (integer *, complex *, integer *, integer *, complex *, complex *, 
			       complex *);
_CPP_PREFIX void CGNLlupqtlsol(integer *, complex *, integer *, integer *, complex *, complex *, 
			       complex *);
_CPP_PREFIX void CGNLlupqusol (integer *, complex *, integer *, integer *, complex *, complex *, 
			       complex *);
_CPP_PREFIX void CGNLlupqtusol(integer *, complex *, integer *, integer *, complex *, complex *, 
			       complex *);
_CPP_PREFIX void CGNLlupqdlsol (integer *, complex *, integer *, integer *, complex *, complex *, 
				complex *);
_CPP_PREFIX void CGNLlupqtdlsol(integer *, complex *, integer *, integer *, complex *, complex *, 
				complex *);
_CPP_PREFIX void CGNLlupqdusol (integer *, complex *, integer *, integer *, complex *, complex *, 
				complex *);
_CPP_PREFIX void CGNLlupqtdusol(integer *, complex *, integer *, integer *, complex *, complex *, 
				complex *);


_CPP_PREFIX void DSPDldlp(integer *, doubleprecision *, integer *, integer *, integer *);
_CPP_PREFIX void DSPDldlpsol (integer *, doubleprecision *, integer *,
			      doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void SSPDldlp(integer *, real *, integer *, integer *, integer *);
_CPP_PREFIX void SSPDldlpsol (integer *, real *, integer *,
			      real *, real *, integer *);

_CPP_PREFIX void ZHPDldlp(integer *, doublecomplex *, integer *, integer *, integer *);
_CPP_PREFIX void ZHPDldlpsol (integer *, doublecomplex *, integer *,
			      doublecomplex *, doublecomplex *, integer *);

_CPP_PREFIX void CHPDldlp(integer *, complex *, integer *, integer *, integer *);
_CPP_PREFIX void CHPDldlpsol (integer *, complex *, integer *,
			      complex *, complex *, integer *);

_CPP_PREFIX integer AMGfactor(SPARSEmat *, AMGlevelmat *, ILUPACKparam *);

_CPP_PREFIX integer DGNLAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX integer DSPDAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX integer DSYMAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);

_CPP_PREFIX integer SGNLAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX integer SSPDAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX integer SSYMAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);

_CPP_PREFIX integer ZGNLAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZHPDAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZHERAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZSYMAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);

_CPP_PREFIX integer CGNLAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CHPDAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CHERAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CSYMAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);



_CPP_PREFIX integer AMGfactor_omp(SPARSEmat *, AMGlevelmat *, ILUPACKparam *);

_CPP_PREFIX integer DGNLAMGfactor_omp(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX integer DSPDAMGfactor_omp(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX integer DSYMAMGfactor_omp(Dmat *, DAMGlevelmat *, DILUPACKparam *);

_CPP_PREFIX integer SGNLAMGfactor_omp(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX integer SSPDAMGfactor_omp(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX integer SSYMAMGfactor_omp(Smat *, SAMGlevelmat *, SILUPACKparam *);

_CPP_PREFIX integer ZGNLAMGfactor_omp(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZHPDAMGfactor_omp(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZHERAMGfactor_omp(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX integer ZSYMAMGfactor_omp(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);

_CPP_PREFIX integer CGNLAMGfactor_omp(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CHPDAMGfactor_omp(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CHERAMGfactor_omp(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX integer CSYMAMGfactor_omp(Cmat *, CAMGlevelmat *, CILUPACKparam *);



_CPP_PREFIX void       DGNLilut(integer *,doubleprecision *,integer *,integer *,integer *,
				doubleprecision *,
				doubleprecision *,integer *,integer *,integer *, 
				doubleprecision *,integer *, integer *);
_CPP_PREFIX void       DGNLilutp(integer *,doubleprecision *,integer *,integer *,integer *,
				 doubleprecision *,doubleprecision *,integer *,
				 doubleprecision *,integer *,integer *,integer *, 
				 doubleprecision *,integer *,integer *,integer *);
_CPP_PREFIX void       DGNLlusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlutsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLludlsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlutdlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLludusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlutdusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlulsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlutlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLluusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void       DGNLlutusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);

_CPP_PREFIX void       SGNLilut(integer *,real *,integer *,integer *,integer *,
				real *,
				real *,integer *,integer *,integer *, 
				real *,integer *, integer *);
_CPP_PREFIX void       SGNLilutp(integer *,real *,integer *,integer *,integer *,
				 real *,real *,integer *,
				 real *,integer *,integer *,integer *, 
				 real *,integer *,integer *,integer *);
_CPP_PREFIX void       SGNLlusol (integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlutsol(integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlulsol (integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlutlsol(integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLluusol (integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlutusol(integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLludlsol (integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlutdlsol(integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLludusol (integer *,real *,real *,real *,integer *,integer *);
_CPP_PREFIX void       SGNLlutdusol(integer *,real *,real *,real *,integer *,integer *);

_CPP_PREFIX void       ZGNLilut (integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,
				 doublecomplex *,integer *,integer *,integer *, 
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLilutp(integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,doubleprecision *,integer *,
				 doublecomplex *,integer *,integer *,integer *, 
				 doublecomplex *,integer *,integer *,integer *);
_CPP_PREFIX void       ZGNLlusol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlutsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlulsol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlutlsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLluusol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlutusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLludlsol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlutdlsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLludusol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void       ZGNLlutdusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,integer *);

_CPP_PREFIX void       CGNLilut (integer *,complex *,integer *,integer *,integer *,
				 real *,
				 complex *,integer *,integer *,integer *, 
				 complex *,integer *,integer *);
_CPP_PREFIX void       CGNLilutp(integer *,complex *,integer *,integer *,integer *,
				 real *,real *,integer *,
				 complex *,integer *,integer *,integer *, 
				 complex *,integer *,integer *,integer *);
_CPP_PREFIX void       CGNLlusol (integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlutsol(integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlulsol (integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlutlsol(integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLluusol (integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlutusol(integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLludlsol (integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlutdlsol(integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLludusol (integer *,complex *,complex *,complex *,integer *,integer *);
_CPP_PREFIX void       CGNLlutdusol(integer *,complex *,complex *,complex *,integer *,integer *);


_CPP_PREFIX void DGNLiluc(integer *,doubleprecision *,integer *,integer *,integer *,
			  doubleprecision *,integer *,doubleprecision *,integer *,integer *,
			  integer *,doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLilucsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			      integer *);
_CPP_PREFIX void DGNLiluctsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			      integer *);
_CPP_PREFIX void DGNLilucdlsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
				integer *);
_CPP_PREFIX void DGNLiluctdlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
				integer *);
_CPP_PREFIX void DGNLilucdusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
				integer *);
_CPP_PREFIX void DGNLiluctdusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
				integer *);
_CPP_PREFIX void DGNLiluclsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			       integer *);
_CPP_PREFIX void DGNLiluctlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			       integer *);
_CPP_PREFIX void DGNLilucusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			       integer *);
_CPP_PREFIX void DGNLiluctusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
			       integer *);

_CPP_PREFIX void DGNLpiluclsol  (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpilucdlsol (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpilucusol  (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpilucdusol (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpiluctlsol (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpiluctdlsol(integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpiluctusol (integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);
_CPP_PREFIX void DGNLpiluctdusol(integer *,integer *, doubleprecision *,doubleprecision *,
				 doubleprecision *,integer *,integer *);

_CPP_PREFIX void DSYMildlc(integer *,doubleprecision *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doubleprecision *,integer *,integer *,
			   doubleprecision *,integer *,integer *,integer *);

_CPP_PREFIX void DSYMildlcnosave(integer *,doubleprecision *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doubleprecision *,integer *,integer *,
				 doubleprecision *,integer *,integer *,integer *,
				 integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSYMildlcsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			      integer *);
_CPP_PREFIX void DSYMpildlcdlsol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSYMpildlcdusol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSYMpildlclsol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSYMpildlcusol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSYMpilucsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			      integer *);
_CPP_PREFIX void DSYMbpilucsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			       integer *,integer *,doubleprecision *);
_CPP_PREFIX void DSYMbpiluclsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			       integer *);
_CPP_PREFIX void DSYMbpilucusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			       integer *);
_CPP_PREFIX void DSSMildlc(integer *,doubleprecision *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doubleprecision *,integer *,integer *,
			   doubleprecision *,integer *,integer *);
_CPP_PREFIX void DSSMildlcnosave(integer *,doubleprecision *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doubleprecision *,integer *,integer *,
				 doubleprecision *,integer *,integer *,
				 integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSSMildlcsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			      integer *);
_CPP_PREFIX void DSSMpildlcdlsol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSSMpildlcdusol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSSMpildlclsol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DSSMpildlcusol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
				 integer *);
_CPP_PREFIX void DGNLpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
			   integer *, doubleprecision *, doubleprecision *, integer*, integer*);
_CPP_PREFIX void DGNLspiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
			    integer *, doubleprecision *, doubleprecision *, integer*, integer*);
_CPP_PREFIX void DGNLmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
			    doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
			    integer *, doubleprecision *, doubleprecision *, integer *, integer*);

_CPP_PREFIX void DSPDpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			   integer *, doubleprecision *, doubleprecision *, 
			   integer *, doubleprecision *,
			   integer *, integer *);

_CPP_PREFIX void DSPDpilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				 doubleprecision *,integer *,integer *,integer *,integer *,
				 doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				 integer *, doubleprecision *, doubleprecision *, 
				 integer *, doubleprecision *,
				 integer *, integer *,
				 integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSPDmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
			    doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			    integer *, doubleprecision *, doubleprecision *, integer *, integer *);

_CPP_PREFIX void DSPDmpilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
				  doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				  integer *, doubleprecision *, doubleprecision *, integer *, integer *,
				  integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSYMpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			   integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
			   doubleprecision *, doubleprecision *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void DSYMpilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				 doubleprecision *,integer *,integer *,integer *,integer *,
				 doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				 integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
				 doubleprecision *, doubleprecision *,
				 integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSYMbpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			    integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
			    doubleprecision *, doubleprecision *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void DSYMbpilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,integer *,integer *,integer *,integer *,
				  doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				  integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
				  doubleprecision *, doubleprecision *, integer *,
				  integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSYMmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			    integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
			    doubleprecision *, doubleprecision *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
/* void DSYMmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doubleprecision *,integer *,integer *,doubleprecision *,integer *,
		integer *, doubleprecision *,  doubleprecision *, integer *, integer *);
*/
_CPP_PREFIX void DSYMmpilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,integer *,integer *,integer *,integer *,
				  doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				  integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
				  doubleprecision *, doubleprecision *,
				  integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void DSYMiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
			  integer *,integer *,integer *,
			  doubleprecision *,integer *,integer *,doubleprecision *,integer *,
			  integer *,integer *);

_CPP_PREFIX void DSYMilucnosave(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
				integer *,integer *,integer *,
				doubleprecision *,integer *,integer *,doubleprecision *,integer *,
				integer *,integer *,
				integer *, doubleprecision *, doubleprecision *, integer *);

_CPP_PREFIX void SGNLiluc(integer *,real *,integer *,integer *,integer *,
			  real *,integer *,real *,integer *,integer *,
			  integer *,real *,integer *,integer *);
_CPP_PREFIX void SGNLilucsol (integer *,real *,real *,real *,integer *,
			      integer *);
_CPP_PREFIX void SGNLiluctsol(integer *,real *,real *,real *,integer *,
			      integer *);
_CPP_PREFIX void SGNLilucdlsol (integer *,real *,real *,real *,integer *,
				integer *);
_CPP_PREFIX void SGNLiluctdlsol(integer *,real *,real *,real *,integer *,
				integer *);
_CPP_PREFIX void SGNLilucdusol (integer *,real *,real *,real *,integer *,
				integer *);
_CPP_PREFIX void SGNLiluctdusol(integer *,real *,real *,real *,integer *,
				integer *);
_CPP_PREFIX void SGNLiluclsol (integer *,real *,real *,real *,integer *,
			       integer *);
_CPP_PREFIX void SGNLiluctlsol(integer *,real *,real *,real *,integer *,
			       integer *);
_CPP_PREFIX void SGNLilucusol (integer *,real *,real *,real *,integer *,
			       integer *);
_CPP_PREFIX void SGNLiluctusol(integer *,real *,real *,real *,integer *,
			       integer *);

_CPP_PREFIX void SGNLpiluclsol  (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpilucdlsol (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpilucusol  (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpilucdusol (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpiluctlsol (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpiluctdlsol(integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpiluctusol (integer *,integer *, real *,real *,
				 real *,integer *,integer *);
_CPP_PREFIX void SGNLpiluctdusol(integer *,integer *, real *,real *,
				 real *,integer *,integer *);

_CPP_PREFIX void SSYMildlc(integer *,real *,integer *,integer *,integer *,
			   real *,integer *,real *,integer *,integer *,
			   real *,integer *,integer *,integer *);
_CPP_PREFIX void SSYMildlcnosave(integer *,real *,integer *,integer *,integer *,
				 real *,integer *,real *,integer *,integer *,
				 real *,integer *,integer *,integer *,
				 integer *, real *, real *, integer *);

_CPP_PREFIX void SSYMildlcsol(integer *,real *,real *,real *,
			      integer *);
_CPP_PREFIX void SSYMpildlcdlsol(integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSYMpildlcdusol(integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSYMpildlclsol (integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSYMpildlcusol (integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSYMpilucsol(integer *,real *,real *,real *,
			      integer *);
_CPP_PREFIX void SSYMbpilucsol(integer *,real *,real *,real *,
			       integer *,integer *,real *);
_CPP_PREFIX void SSYMbpiluclsol(integer *,real *,real *,real *, integer *);
_CPP_PREFIX void SSYMbpilucusol(integer *,real *,real *,real *, integer *);

_CPP_PREFIX void SSSMildlc(integer *,real *,integer *,integer *,integer *,
			   real *,integer *,real *,integer *,integer *,
			   real *,integer *,integer *);

_CPP_PREFIX void SSSMildlcnosave(integer *,real *,integer *,integer *,integer *,
				 real *,integer *,real *,integer *,integer *,
				 real *,integer *,integer *,
				 integer *, real *, real *, integer *);

_CPP_PREFIX void SSSMildlcsol(integer *,real *,real *,real *,
			      integer *);
_CPP_PREFIX void SSSMpildlcdlsol(integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSSMpildlcdusol(integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSSMpildlclsol (integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SSSMpildlcusol (integer *,integer *,real *,real *,real *,
				 integer *);
_CPP_PREFIX void SGNLpiluc(integer *,real *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   real *,integer *,integer *,integer *,real *,integer *,
			   integer *, real *, real *, integer *, integer*);
_CPP_PREFIX void SGNLspiluc(integer *,real *,integer *,integer *,integer *,real *,
			    real *,integer *,integer *,integer *,integer *,
			    real *,integer *,integer *,integer *,real *,integer *,
			    integer *, real *, real *, integer *, integer*);
_CPP_PREFIX void SGNLmpiluc(integer *,real *,integer *,integer *,integer *,real *,
			    real *,real *,integer *,integer *,integer *,integer *,
			    real *,integer *,integer *,integer *,real *,integer *,
			    integer *, real *, real *, integer *, integer*);

_CPP_PREFIX void SSPDpiluc(integer *,real *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   real *,integer *,integer *,real *,integer *,
			   integer *, real *, real *, 
			   integer *, real *,
			   integer *, integer *);

_CPP_PREFIX void SSPDpilucnosave(integer *,real *,integer *,integer *,integer *,real *,
				 real *,integer *,integer *,integer *,integer *,
				 real *,integer *,integer *,real *,integer *,
				 integer *, real *, real *, 
				 integer *, real *,
				 integer *, integer *,
				 integer *, real *, real *, integer *);

_CPP_PREFIX void SSPDmpiluc(integer *,real *,integer *,integer *,integer *,real *,
			    real *,real *,integer *,integer *,integer *,integer *,
			    real *,integer *,integer *,real *,integer *,
			    integer *, real *, real *, integer *, integer *);

_CPP_PREFIX void SSPDmpilucnosave(integer *,real *,integer *,integer *,integer *,real *,
				  real *,real *,integer *,integer *,integer *,integer *,
				  real *,integer *,integer *,real *,integer *,
				  integer *, real *, real *, integer *, integer *,
				  integer *, real *, real *, integer *);

_CPP_PREFIX void SSYMpiluc(integer *,real *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   real *,integer *,integer *,real *,integer *,
			   integer *, real *, real *, integer *, integer *, integer *,
			   real *, real *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void SSYMpilucnosave(integer *,real *,integer *,integer *,integer *,real *,
				 real *,integer *,integer *,integer *,integer *,
				 real *,integer *,integer *,real *,integer *,
				 integer *, real *, real *, integer *, integer *, integer *,
				 real *, real *,
				 integer *, real *, real *, integer *);

_CPP_PREFIX void SSYMbpiluc(integer *,real *,integer *,integer *,integer *,real *,
			    real *,integer *,integer *,integer *,integer *,
			    real *,integer *,integer *,real *,integer *,
			    integer *, real *, real *, integer *, integer *, integer *,
			    real *, real *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
/* void SSYMmpiluc(integer *,real *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		real *,integer *,integer *,real *,integer *,
		integer *, real *, real *, integer *, integer *);
*/
_CPP_PREFIX void SSYMbpilucnosave(integer *,real *,integer *,integer *,integer *,real *,
				  real *,integer *,integer *,integer *,integer *,
				  real *,integer *,integer *,real *,integer *,
				  integer *, real *, real *, integer *, integer *, integer *,
				  real *, real *, integer *,
				  integer *, real *, real *, integer *);

_CPP_PREFIX void SSYMiluc(integer *,real *,integer *,integer *,integer *,real *,
			  integer *,integer *,integer *,
			  real *,integer *,integer *,real *,integer *,
			  integer *, integer *);

_CPP_PREFIX void SSYMilucnosave(integer *,real *,integer *,integer *,integer *,real *,
				integer *,integer *,integer *,
				real *,integer *,integer *,real *,integer *,
				integer *, integer *,
				integer *, real *, real *, integer *);

_CPP_PREFIX void ZGNLiluc(integer *,doublecomplex *,integer *,integer *,integer *,
			  doubleprecision *,integer *,doublecomplex *,integer *,integer *,
			  integer *,doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLilucsol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			      integer *);
_CPP_PREFIX void ZGNLiluctsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			      integer *);
_CPP_PREFIX void ZGNLilucdlsol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
				integer *);
_CPP_PREFIX void ZGNLiluctdlsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
				integer *);
_CPP_PREFIX void ZGNLilucdusol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
				integer *);
_CPP_PREFIX void ZGNLiluctdusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
				integer *);
_CPP_PREFIX void ZGNLiluclsol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			       integer *);
_CPP_PREFIX void ZGNLiluctlsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			       integer *);
_CPP_PREFIX void ZGNLilucusol (integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			       integer *);
_CPP_PREFIX void ZGNLiluctusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,integer *,
			       integer *);

_CPP_PREFIX void ZGNLpiluclsol  (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpilucdlsol (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpilucusol  (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpilucdusol (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpiluctlsol (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpiluctdlsol(integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpiluctusol (integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);
_CPP_PREFIX void ZGNLpiluctdusol(integer *,integer *, doublecomplex *,doublecomplex *,
				 doublecomplex *,integer *,integer *);

_CPP_PREFIX void ZHERildlc(integer *,doublecomplex *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doublecomplex *,integer *,integer *,
			   doublecomplex *,integer *,integer *,integer *);

_CPP_PREFIX void ZHERildlcnosave(integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doublecomplex *,integer *,integer *,
				 doublecomplex *,integer *,integer *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *);


_CPP_PREFIX void ZHERildlcsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZHERpildlcdlsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZHERpildlcdusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZHERpildlclsol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZHERpildlcusol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZHERpilucsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZHERbpilucsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			       integer *,integer *,doublecomplex *);
_CPP_PREFIX void ZHERbpiluclsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZHERbpilucusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZSYMpilucsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZSYMbpilucsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			       integer *, integer *,doublecomplex *);
_CPP_PREFIX void ZSYMbpiluclsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZSYMbpilucusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);

_CPP_PREFIX void ZSYMildlc(integer *,doublecomplex *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doublecomplex *,integer *,integer *,
			   doublecomplex *,integer *,integer *,integer *);

_CPP_PREFIX void ZSYMildlcnosave(integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doublecomplex *,integer *,integer *,
				 doublecomplex *,integer *,integer *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZSYMildlcsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZSYMpildlcdlsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSYMpildlcdusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSYMpildlclsol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSYMpildlcusol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);

_CPP_PREFIX void ZSHRildlc(integer *,doublecomplex *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doublecomplex *,integer *,integer *,
			   doublecomplex *,integer *,integer *);

_CPP_PREFIX void ZSHRildlcnosave(integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doublecomplex *,integer *,integer *,
				 doublecomplex *,integer *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZSHRildlcsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZSHRpildlcdlsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSHRpildlcdusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSHRpildlclsol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSHRpildlcusol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);

_CPP_PREFIX void ZSSMildlc(integer *,doublecomplex *,integer *,integer *,integer *,
			   doubleprecision *,integer *,doublecomplex *,integer *,integer *,
			   doublecomplex *,integer *,integer *);

_CPP_PREFIX void ZSSMildlcnosave(integer *,doublecomplex *,integer *,integer *,integer *,
				 doubleprecision *,integer *,doublecomplex *,integer *,integer *,
				 doublecomplex *,integer *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZSSMildlcsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *,
			      integer *);
_CPP_PREFIX void ZSSMpildlcdlsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSSMpildlcdusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSSMpildlclsol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZSSMpildlcusol (integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *,
				 integer *);
_CPP_PREFIX void ZGNLpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doublecomplex *,integer *,integer *,integer *,doublecomplex *,integer *,
			   integer *, doubleprecision *, doublecomplex *, integer *, integer *);
_CPP_PREFIX void ZGNLspiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doublecomplex *,integer *,integer *,integer *,doublecomplex *,integer *,
			    integer *, doubleprecision *, doublecomplex *, integer *, integer *);
_CPP_PREFIX void ZGNLmpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
			    doublecomplex *,integer *,integer *,integer *,doublecomplex *,integer *,
			    integer *, doubleprecision *, doublecomplex *, integer *, integer *);

_CPP_PREFIX void ZHPDpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			   integer *, doubleprecision *, doublecomplex *,
			   integer *, doublecomplex *, 
			   integer *, integer *);

_CPP_PREFIX void ZHPDpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				 doubleprecision *,integer *,integer *,integer *,integer *,
				 doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				 integer *, doubleprecision *, doublecomplex *,
				 integer *, doublecomplex *, 
				 integer *, integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZHPDmpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
			    doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			    integer *, doubleprecision *, doublecomplex *, integer *, integer *);

_CPP_PREFIX void ZHPDmpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
				  doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *, integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZHERpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			   integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
			   doubleprecision *, doubleprecision *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void ZHERpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				 doubleprecision *,integer *,integer *,integer *,integer *,
				 doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
				 doubleprecision *, doubleprecision *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZHERbpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			    integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
			    doubleprecision *, doubleprecision *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
/* void ZHERmpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doublecomplex *,integer *,integer *,doublecomplex *,integer *,
		integer *, doubleprecision *, doublecomplex *, integer *, integer *);
*/
_CPP_PREFIX void ZHERbpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,integer *,integer *,integer *,integer *,
				  doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
				  doubleprecision *, doubleprecision *, integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZHERiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			  integer *,integer *,integer *,
			  doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			  integer *,integer *);

_CPP_PREFIX void ZHERilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				integer *,integer *,integer *,
				doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				integer *,integer *,
				integer *, doubleprecision *, doublecomplex *, integer *);


_CPP_PREFIX void ZSYMpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			   doubleprecision *,integer *,integer *,integer *,integer *,
			   doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			   integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
			   doubleprecision *, doubleprecision *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void ZSYMpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				 doubleprecision *,integer *,integer *,integer *,integer *,
				 doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				 integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
				 doubleprecision *, doubleprecision *,
				 integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZSYMbpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			    doubleprecision *,integer *,integer *,integer *,integer *,
			    doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			    integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
			    doubleprecision *, doubleprecision *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
/*void ZSYMmpiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doublecomplex *,integer *,integer *,doublecomplex *,integer *,
		integer *, doubleprecision *, doublecomplex *, integer *, integer *);
*/
_CPP_PREFIX void ZSYMbpilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				  doubleprecision *,integer *,integer *,integer *,integer *,
				  doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *, integer *, integer *,
				  doubleprecision *, doubleprecision *, integer *,
				  integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void ZSYMiluc(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
			  integer *,integer *,integer *,
			  doublecomplex *,integer *,integer *,doublecomplex *,integer *,
			  integer *, integer *);

_CPP_PREFIX void ZSYMilucnosave(integer *,doublecomplex *,integer *,integer *,integer *,doubleprecision *,
				integer *,integer *,integer *,
				doublecomplex *,integer *,integer *,doublecomplex *,integer *,
				integer *, integer *,
				integer *, doubleprecision *, doublecomplex *, integer *);

_CPP_PREFIX void CGNLiluc(integer *,complex *,integer *,integer *,integer *,
			  real *,integer *,complex *,integer *,integer *,
			  integer *,complex *,integer *,integer *);
_CPP_PREFIX void CGNLilucsol (integer *,complex *,complex *,complex *,integer *,
			      integer *);
_CPP_PREFIX void CGNLiluctsol(integer *,complex *,complex *,complex *,integer *,
			      integer *);
_CPP_PREFIX void CGNLilucdlsol (integer *,complex *,complex *,complex *,integer *,
				integer *);
_CPP_PREFIX void CGNLiluctdlsol(integer *,complex *,complex *,complex *,integer *,
				integer *);
_CPP_PREFIX void CGNLilucdusol (integer *,complex *,complex *,complex *,integer *,
				integer *);
_CPP_PREFIX void CGNLiluctdusol(integer *,complex *,complex *,complex *,integer *,
				integer *);
_CPP_PREFIX void CGNLilucusol (integer *,complex *,complex *,complex *,integer *,
			       integer *);
_CPP_PREFIX void CGNLiluctusol(integer *,complex *,complex *,complex *,integer *,
			       integer *);
_CPP_PREFIX void CGNLiluclsol (integer *,complex *,complex *,complex *,integer *,
			       integer *);
_CPP_PREFIX void CGNLiluctlsol(integer *,complex *,complex *,complex *,integer *,
			       integer *);

_CPP_PREFIX void CGNLpiluclsol  (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpilucdlsol (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpilucusol  (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpilucdusol (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpiluctlsol (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpiluctdlsol(integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpiluctusol (integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);
_CPP_PREFIX void CGNLpiluctdusol(integer *,integer *, complex *,complex *,
				 complex *,integer *,integer *);

_CPP_PREFIX void CHERildlc(integer *,complex *,integer *,integer *,integer *,
			   real *,integer *,complex *,integer *,integer *,
			   complex *,integer *,integer *,integer *);

_CPP_PREFIX void CHERildlcnosave(integer *,complex *,integer *,integer *,integer *,
				 real *,integer *,complex *,integer *,integer *,
				 complex *,integer *,integer *,integer *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CHERildlcsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CHERpildlcdlsol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CHERpildlcdusol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CHERpildlclsol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CHERpildlcusol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CHERpilucsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CSYMpilucsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CHERbpilucsol(integer *,complex *,complex *,complex *,
			       integer *,integer *,complex *);
_CPP_PREFIX void CHERbpiluclsol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERbpilucusol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMbpilucsol(integer *,complex *,complex *,complex *,
			       integer *, integer *,complex *);
_CPP_PREFIX void CSYMbpiluclsol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMbpilucusol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   complex *,integer *,integer *,complex *,integer *,
			   integer *, real *, complex *, integer *, integer *, integer *,
			   real *, real *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void CSYMpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				 real *,integer *,integer *,integer *,integer *,
				 complex *,integer *,integer *,complex *,integer *,
				 integer *, real *, complex *, integer *, integer *, integer *,
				 real *, real *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CSYMbpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			    real *,integer *,integer *,integer *,integer *,
			    complex *,integer *,integer *,complex *,integer *,
			    integer *, real *, complex *, integer *, integer *, integer *,
			    real *, real *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void CSYMbpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				  real *,integer *,integer *,integer *,integer *,
				  complex *,integer *,integer *,complex *,integer *,
				  integer *, real *, complex *, integer *, integer *, integer *,
				  real *, real *, integer *,
				  integer *, real *, complex *, integer *);

_CPP_PREFIX void CHERpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   complex *,integer *,integer *,complex *,integer *,
			   integer *, real *, complex *, integer *, integer *, integer *,
			   real *, real *);
/*
			   integer *,integer *,integer *,integer *,integer *,integer *,
			   integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void CHERpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				 real *,integer *,integer *,integer *,integer *,
				 complex *,integer *,integer *,complex *,integer *,
				 integer *, real *, complex *, integer *, integer *, integer *,
				 real *, real *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CHERbpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			    real *,integer *,integer *,integer *,integer *,
			    complex *,integer *,integer *,complex *,integer *,
			    integer *, real *, complex *, integer *, integer *, integer *,
			    real *, real *, integer *);
/*
			    integer *,integer *,integer *,integer *,integer *,integer *,
			    integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
*/
_CPP_PREFIX void CHERbpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				  real *,integer *,integer *,integer *,integer *,
				  complex *,integer *,integer *,complex *,integer *,
				  integer *, real *, complex *, integer *, integer *, integer *,
				  real *, real *, integer *,
				  integer *, real *, complex *, integer *);

_CPP_PREFIX void CSYMiluc(integer *,complex *,integer *,integer *,integer *,real *,
			  integer *,integer *,integer *,
			  complex *,integer *,integer *,complex *,integer *,
			  integer *, integer *);

_CPP_PREFIX void CSYMilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				integer *,integer *,integer *,
				complex *,integer *,integer *,complex *,integer *,
				integer *, integer *,
				integer *, real *, complex *, integer *);

_CPP_PREFIX void CHERiluc(integer *,complex *,integer *,integer *,integer *,real *,
			  integer *,integer *,integer *,
			  complex *,integer *,integer *,complex *,integer *,
			  integer *, integer *);

_CPP_PREFIX void CHERilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				integer *,integer *,integer *,
				complex *,integer *,integer *,complex *,integer *,
				integer *, integer *,
				integer *, real *, complex *, integer *);

_CPP_PREFIX void CSYMildlc(integer *,complex *,integer *,integer *,integer *,
			   real *,integer *,complex *,integer *,integer *,
			   complex *,integer *,integer *,integer *);

_CPP_PREFIX void CSYMildlcnosave(integer *,complex *,integer *,integer *,integer *,
				 real *,integer *,complex *,integer *,integer *,
				 complex *,integer *,integer *,integer *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CSYMildlcsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CSYMpildlcdlsol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSYMpildlcdusol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSYMpildlclsol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSYMpildlcusol (integer *,integer *,complex *,complex *,complex *,
				 integer *);

_CPP_PREFIX void CSHRildlc(integer *,complex *,integer *,integer *,integer *,
			   real *,integer *,complex *,integer *,integer *,
			   complex *,integer *,integer *);

_CPP_PREFIX void CSHRildlcnosave(integer *,complex *,integer *,integer *,integer *,
				 real *,integer *,complex *,integer *,integer *,
				 complex *,integer *,integer *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CSHRildlcsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CSHRpildlcdlsol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSHRpildlcdusol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSHRpildlclsol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSHRpildlcusol (integer *,integer *,complex *,complex *,complex *,
				 integer *);

_CPP_PREFIX void CSSMildlc(integer *,complex *,integer *,integer *,integer *,
			   real *,integer *,complex *,integer *,integer *,
			   complex *,integer *,integer *);

_CPP_PREFIX void CSSMildlcnosave(integer *,complex *,integer *,integer *,integer *,
				 real *,integer *,complex *,integer *,integer *,
				 complex *,integer *,integer *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CSSMildlcsol(integer *,complex *,complex *,complex *,
			      integer *);
_CPP_PREFIX void CSSMpildlcdlsol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSSMpildlcdusol(integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSSMpildlclsol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CSSMpildlclsol (integer *,integer *,complex *,complex *,complex *,
				 integer *);
_CPP_PREFIX void CGNLpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   complex *,integer *,integer *,integer *,complex *,integer *,
			   integer *, real *, complex *, integer *, integer *);
_CPP_PREFIX void CGNLspiluc(integer *,complex *,integer *,integer *,integer *,real *,
			    real *,integer *,integer *,integer *,integer *,
			    complex *,integer *,integer *,integer *,complex *,integer *,
			    integer *, real *, complex *, integer *, integer *);
_CPP_PREFIX void CGNLmpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			    real *,real *,integer *,integer *,integer *,integer *,
			    complex *,integer *,integer *,integer *,complex *,integer *,
			    integer *, real *, complex *, integer *, integer *);

_CPP_PREFIX void CHPDpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			   real *,integer *,integer *,integer *,integer *,
			   complex *,integer *,integer *,complex *,integer *,
			   integer *, real *, complex *, 
			   integer *, complex *, 
			   integer *, integer *);

_CPP_PREFIX void CHPDpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				 real *,integer *,integer *,integer *,integer *,
				 complex *,integer *,integer *,complex *,integer *,
				 integer *, real *, complex *, 
				 integer *, complex *, 
				 integer *, integer *,
				 integer *, real *, complex *, integer *);

_CPP_PREFIX void CHPDmpiluc(integer *,complex *,integer *,integer *,integer *,real *,
			    real *,real *,integer *,integer *,integer *,integer *,
			    complex *,integer *,integer *,complex *,integer *,
			    integer *, real *, complex *, integer *, integer *);

_CPP_PREFIX void CHPDmpilucnosave(integer *,complex *,integer *,integer *,integer *,real *,
				  real *,real *,integer *,integer *,integer *,integer *,
				  complex *,integer *,integer *,complex *,integer *,
				  integer *, real *, complex *, integer *, integer *,
				  integer *, real *, complex *, integer *);


/* *********************************************** */
/* ******      Definitions for orderings     ***** */

_CPP_PREFIX integer    DGNLperm_null        (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cnull        (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_nd          (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_nd_fc       (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_nd_fcv      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_nd_fc       (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_nd_fcv      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_rcm         (Dmat, doubleprecision *,doubleprecision *,
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_rcm_fc      (Dmat, doubleprecision *,doubleprecision *,
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_rcm_fcv     (Dmat, doubleprecision *,doubleprecision *,
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mmd         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mmd_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mmd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mmd_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mmd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amf         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amf_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_amf_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_amf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camf         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camf_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_camf_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_camf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amd         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camd        (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amd_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camd_fc     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_amd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_camd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_amd_fc      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_camd_fc     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_amd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_camd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_metis_e     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mtmetis     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_metis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_metis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_metis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_metis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_e     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_cmetis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_cmetis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_metis_n     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);

_CPP_PREFIX integer    DGNLpartition_metis_n(Dmat, integer *, integer *, integer, integer,
					     integer *, integer, DILUPACKparam *);


_CPP_PREFIX integer    DGNLperm_metis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_metis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_cmetis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_metis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mtmetis     (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_cmetis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_metis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_cmetis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_pq          (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_fc          (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_fc          (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_fcv         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_fcv         (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_p           (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_indset      (Dmat, doubleprecision *,doubleprecision *, 
					     integer *,integer *, integer *, DILUPACKparam *);


_CPP_PREFIX integer    DGNLperm_mwm_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camd       (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwm_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);

_CPP_PREFIX integer    DSYMperm_mwm_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_rcm_sp     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_mmd_sp     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amf_sp     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camf_sp     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camd       (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amd_sp     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camd_sp    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mwm_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);


_CPP_PREFIX integer    DGNLperm_matching_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amf        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camf        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amd        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camd       (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_matching_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);

_CPP_PREFIX integer    DSYMperm_matching_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_rcm_sp     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_mmd_sp        (Dmat, doubleprecision *,doubleprecision *, 
							integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amf        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amf_sp     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camf        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camf_sp     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amd        (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camd       (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amd_sp     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camd_sp    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_matching_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						     integer *,integer *, integer *, DILUPACKparam *);


_CPP_PREFIX integer    DSYMperm_mc64_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_mtmetis        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_rcm_sp        (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_mmd_sp        (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amf        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amf_sp        (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camf        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camf_sp        (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amd        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camd       (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amd_sp     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camd_sp    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_n_sp    (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_n_sp    (Dmat, doubleprecision *,doubleprecision *, 
						    integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DSYMperm_mc64_null        (Dmat, doubleprecision *,doubleprecision *, 
						  integer *,integer *, integer *, DILUPACKparam *);



_CPP_PREFIX integer    DGNLperm_mc64_null       (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_mtmetis       (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_amf        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_camf        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_amd        (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_camd       (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mc64_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						 integer *,integer *, integer *, DILUPACKparam *);

_CPP_PREFIX integer    DGNLperm_mwa_rcm        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_mmd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camf        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camf_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amd        (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camd       (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camd_fc    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_camd_fcv   (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_e    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_n    (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);
_CPP_PREFIX integer    DGNLperm_mwa_cmetis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
						integer *,integer *, integer *, DILUPACKparam *);


_CPP_PREFIX integer    SGNLperm_null        (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cnull       (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_nd          (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_nd_fc       (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_nd_fcv      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_nd_fc       (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_nd_fcv      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_rcm         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_rcm_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_rcm_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mmd         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mmd_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mmd_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mmd_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mmd_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amf         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amf_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amf_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_amf_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_amf_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camf         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camf_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camf_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_camf_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_camf_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amd         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camd        (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amd_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camd_fc     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_amd_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_camd_fcv    (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_amd_fc      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_camd_fc     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_amd_fcv     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_camd_fcv    (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_metis_e     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mtmetis     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_metis_e_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_metis_e_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_metis_e_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_metis_e_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_e     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_e_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_e_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_cmetis_e_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_cmetis_e_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_metis_n     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_n     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);

_CPP_PREFIX integer    SGNLpartition_metis_n(Smat, integer *, integer *, integer, integer,
					     integer *, integer, SILUPACKparam *);



_CPP_PREFIX integer    SGNLperm_metis_n_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_n_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_metis_n_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_cmetis_n_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_metis_n_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mtmetis     (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_cmetis_n_fc  (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_metis_n_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_cmetis_n_fcv (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_pq          (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_fc          (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_fc          (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_fcv         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_fcv         (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_indset      (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_p           (Smat, real *,real *, integer *,integer *,
					     integer *, SILUPACKparam *);

_CPP_PREFIX integer    SGNLperm_mwm_rcm        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_rcm_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_rcm_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_mmd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_mmd_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_mmd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amf_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camf_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camd       (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amd_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camd_sp    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camd_fc    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_amd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_camd_fcv   (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_e_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_e_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_n_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_n_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwm_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);

_CPP_PREFIX integer    SSYMperm_mwm_rcm        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_rcm_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_rcm_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_mmd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_mmd_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_mmd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amf_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camf_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camd       (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amd_sp     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camd_sp    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camd_fc    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_amd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_camd_fcv   (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_e_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_e_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_n_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_n_sp (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mwm_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);


_CPP_PREFIX integer    SGNLperm_matching_rcm        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_rcm_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_mmd        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_mmd_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amf        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amf_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amf_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camf        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camf_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camf_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amd        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camd       (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amd_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camd_fc    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_amd_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_camd_fcv   (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_e    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_e_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_e    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_n    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_n    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_n_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_matching_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);

_CPP_PREFIX integer    SSYMperm_matching_rcm        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_rcm_sp     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_rcm_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_mmd        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_mmd_sp     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_mmd_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amf        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amf_sp     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amf_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amf_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camf        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camf_sp     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camf_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camf_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amd        (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camd       (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amd_sp     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camd_sp    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amd_fc     (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camd_fc    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_amd_fcv    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_camd_fcv   (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_e    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_e_sp (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_e_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_e    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_e_sp (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_n    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_n    (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_n_sp (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_n_sp (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_n_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_matching_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						     integer *, SILUPACKparam *);


_CPP_PREFIX integer    SSYMperm_mc64_rcm        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_mtmetis        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_rcm_sp        (Smat, real *,real *, integer *,integer *,
						    integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_rcm_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_mmd        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_mmd_sp        (Smat, real *,real *, integer *,integer *,
						    integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_mmd_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amf        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amf_sp     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amf_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amf_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camf        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camf_sp     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camf_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camf_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amd        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camd       (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amd_sp     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camd_sp    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amd_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camd_fc    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_amd_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_camd_fcv   (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_e    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_e_sp (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_e_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_e    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_e_sp (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_n    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_n    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_n_sp (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_n_sp (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_n_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSYMperm_mc64_null        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);



_CPP_PREFIX integer    SGNLperm_mc64_null       (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_mtmetis    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_rcm        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_rcm_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_mmd        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_mmd_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amf        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amf_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amf_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camf        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camf_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camf_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amd        (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camd       (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amd_fc     (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camd_fc    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_amd_fcv    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_camd_fcv   (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_e    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_e_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_e    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_n    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_n    (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_n_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mc64_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						 integer *, SILUPACKparam *);

_CPP_PREFIX integer    SGNLperm_mwa_rcm        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_rcm_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_rcm_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_mmd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_mmd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_mmd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camf        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camf_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camf_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amd        (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camd       (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amd_fc     (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camd_fc    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_amd_fcv    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_camd_fcv   (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_e    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_e_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_e_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_n    (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_n_fc (Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_metis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);
_CPP_PREFIX integer    SGNLperm_mwa_cmetis_n_fcv(Smat, real *,real *, integer *,integer *,
						integer *, SILUPACKparam *);


_CPP_PREFIX integer    ZGNLperm_null        (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cnull        (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_nd          (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_nd_fc       (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_nd_fcv      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_nd_fc       (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_nd_fcv      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_rcm         (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_rcm_fc      (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_rcm_fcv     (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mmd         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mmd_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mmd_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mmd_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mmd_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amf         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amf_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amf_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_amf_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_amf_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camf         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camf_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camf_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_camf_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_camf_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amd         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camd        (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amd_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camd_fc     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_amd_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_camd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_amd_fc      (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_camd_fc     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_amd_fcv     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_camd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_metis_e     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mtmetis     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_metis_e_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_metis_e_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_metis_e_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_metis_e_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_e     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_e_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_e_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_cmetis_e_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_cmetis_e_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_metis_n     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_n     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZGNLpartition_metis_n(Zmat, integer *, integer *, integer, integer,
					     integer *, integer, ZILUPACKparam *);


_CPP_PREFIX integer    ZGNLperm_metis_n_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_n_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_metis_n_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_cmetis_n_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_metis_n_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mtmetis     (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_cmetis_n_fc  (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_metis_n_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_cmetis_n_fcv (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_pq          (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_fc          (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_fc          (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_fcv         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_fcv         (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_indset      (Zmat, doublecomplex *,doublecomplex *,
					     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_p           (Zmat, doublecomplex *,doublecomplex *, 
					     integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZGNLperm_mwm_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amd        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camd       (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwm_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZHERperm_mtmetis        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_mmd_sp     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amf_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camf_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amd        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camd       (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mwm_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZSYMperm_mwm_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_mmd_sp     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amf_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camf_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amd        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camd       (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mwm_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);


_CPP_PREFIX integer    ZGNLperm_matching_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amd        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camd       (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_matching_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZHERperm_matching_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_mmd_sp     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amf_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camf_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amd        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camd       (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_matching_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZSYMperm_matching_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_mmd_sp     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amf_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camf        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camf_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amd        (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camd       (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_matching_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						     integer *,integer *, integer *, ZILUPACKparam *);


_CPP_PREFIX integer    ZHERperm_mc64_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_mtmetis    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_mmd_sp        (Zmat, doublecomplex *,doublecomplex *, 
						    integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amf_sp        (Zmat, doublecomplex *,doublecomplex *,
						    integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camf_sp        (Zmat, doublecomplex *,doublecomplex *,
						    integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amd        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camd       (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_e_sp    (Zmat, doublecomplex *,doublecomplex *,
						    integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_e_sp    (Zmat, doublecomplex *,doublecomplex *,
						    integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZHERperm_mc64_null        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZSYMperm_mc64_rcm        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_mtmetis    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_rcm_sp     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_rcm_fc     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_mmd        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_mmd_sp     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_mmd_fc     (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amf_sp     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camf_sp     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amd        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camd       (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amd_sp     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camd_sp    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_e_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_n_sp (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZSYMperm_mc64_null        (Zmat, doublecomplex *,doublecomplex *, 
						 integer *,integer *, integer *, ZILUPACKparam *);


_CPP_PREFIX integer    ZGNLperm_mc64_null       (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_mtmetis    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cnull       (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_rcm        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_rcm_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_mmd        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_mmd_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camf        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amd        (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camd       (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mc64_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						 integer *,integer *, integer *, ZILUPACKparam *);


_CPP_PREFIX integer    ZGNLperm_mwa_rcm        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_rcm_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_rcm_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_mmd        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_mmd_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_mmd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camf        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camf_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camf_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amd        (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camd       (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amd_fc     (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camd_fc    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_amd_fcv    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_camd_fcv   (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_e    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_e_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_e_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_n    (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_n_fc (Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_metis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);
_CPP_PREFIX integer    ZGNLperm_mwa_cmetis_n_fcv(Zmat, doublecomplex *,doublecomplex *,
						integer *,integer *, integer *, ZILUPACKparam *);


_CPP_PREFIX integer    CGNLperm_null        (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cnull        (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_nd          (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_nd_fc       (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_nd_fcv      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_nd_fc       (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_nd_fcv      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_rcm         (Cmat, complex *,complex *, integer *,integer *, 
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_rcm_fc      (Cmat, complex *,complex *, integer *,integer *, 
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_rcm_fcv     (Cmat, complex *,complex *, integer *,integer *, 
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mmd         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mmd_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mmd_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mmd_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mmd_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amf         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amf_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amf_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_amf_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_amf_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camf         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camf_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camf_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_camf_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_camf_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amd         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camd        (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amd_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camd_fc     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_amd_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_camd_fcv    (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_amd_fc      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_camd_fc     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_amd_fcv     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_camd_fcv    (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_metis_e     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mtmetis     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_metis_e_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_metis_e_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_metis_e_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mtmetis     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_metis_e_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_e     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_e_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_e_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_cmetis_e_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_cmetis_e_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_metis_n     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_n     (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);

_CPP_PREFIX integer    CGNLpartition_metis_n(Cmat, integer *, integer *, integer, integer,
					     integer *, integer, CILUPACKparam *);

_CPP_PREFIX integer    CGNLperm_metis_n_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_n_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_metis_n_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_cmetis_n_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_metis_n_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_cmetis_n_fc  (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_metis_n_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_cmetis_n_fcv (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_pq          (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_fc          (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_fc          (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_fcv         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_fcv         (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_indset      (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_p           (Cmat, complex *,complex *, integer *,integer *,
					     integer *, CILUPACKparam *);

_CPP_PREFIX integer    CGNLperm_mwm_rcm        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_mmd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camd       (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwm_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);

_CPP_PREFIX integer    CHERperm_mtmetis        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_rcm        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_mmd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amf_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camf_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camd       (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mwm_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);

_CPP_PREFIX integer    CSYMperm_mwm_rcm        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_mmd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amf_sp    (Cmat, complex *,complex *, integer *,integer *,
					       integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camf_sp    (Cmat, complex *,complex *, integer *,integer *,
					       integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camd       (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mwm_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);

_CPP_PREFIX integer    CGNLperm_matching_rcm        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_mmd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camd       (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_matching_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);

_CPP_PREFIX integer    CHERperm_matching_rcm        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_mmd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amf_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camf_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camd       (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_matching_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);

_CPP_PREFIX integer    CSYMperm_matching_rcm        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_mmd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amf_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camf        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camf_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amd        (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camd       (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_matching_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						     integer *, CILUPACKparam *);


_CPP_PREFIX integer    CHERperm_mc64_rcm        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_mtmetis    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_mmd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amf_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camf_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camd       (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CHERperm_mc64_null        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);

_CPP_PREFIX integer    CSYMperm_mc64_rcm        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_mtmetis    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_rcm_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_mmd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_mmd_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amf_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camf_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camd       (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amd_sp     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camd_sp    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_e_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_n_sp (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CSYMperm_mc64_null        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);


_CPP_PREFIX integer    CGNLperm_mc64_null       (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_mtmetis    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_rcm        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_mmd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camf        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amd        (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camd       (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mc64_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						 integer *, CILUPACKparam *);


_CPP_PREFIX integer    CGNLperm_mwa_rcm        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_rcm_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_rcm_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_mmd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_mmd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_mmd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camf        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camf_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camf_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amd        (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camd       (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amd_fc     (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camd_fc    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_amd_fcv    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_camd_fcv   (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_e    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_e_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_e_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_n    (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_n_fc (Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_metis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);
_CPP_PREFIX integer    CGNLperm_mwa_cmetis_n_fcv(Cmat, complex *,complex *, integer *,integer *,
						integer *, CILUPACKparam *);




#define DSPDperm_null  DGNLperm_null
#define SSPDperm_null  SGNLperm_null
#define ZHPDperm_null  ZGNLperm_null
#define CHPDperm_null  CGNLperm_null

#define DSPDpermnull  DGNLperm_null
#define SSPDpermnull  SGNLperm_null
#define ZHPDpermnull  ZGNLperm_null
#define CHPDpermnull  CGNLperm_null

#define DSYMperm_null  DGNLperm_null
#define SSYMperm_null  SGNLperm_null
#define ZSYMperm_null  ZGNLperm_null
#define CSYMperm_null  CGNLperm_null
#define ZHERperm_null  ZGNLperm_null
#define CHERperm_null  CGNLperm_null

#define DSPDperm_cnull  DGNLperm_cnull
#define SSPDperm_cnull  SGNLperm_cnull
#define ZHPDperm_cnull  ZGNLperm_cnull
#define CHPDperm_cnull  CGNLperm_cnull

#define DSPDpermcnull  DGNLperm_cnull
#define SSPDpermcnull  SGNLperm_cnull
#define ZHPDpermcnull  ZGNLperm_cnull
#define CHPDpermcnull  CGNLperm_cnull

#define DSYMperm_cnull  DGNLperm_cnull
#define SSYMperm_cnull  SGNLperm_cnull
#define ZSYMperm_cnull  ZGNLperm_cnull
#define CSYMperm_cnull  CGNLperm_cnull
#define ZHERperm_cnull  ZGNLperm_cnull
#define CHERperm_cnull  CGNLperm_cnull


#define DSPDperm_nd    DGNLperm_nd
#define SSPDperm_nd    SGNLperm_nd
#define ZHPDperm_nd    ZGNLperm_nd
#define CHPDperm_nd    CGNLperm_nd

#define DSPDpermnd    DGNLperm_nd
#define SSPDpermnd    SGNLperm_nd
#define ZHPDpermnd    ZGNLperm_nd
#define CHPDpermnd    CGNLperm_nd

#define DSYMperm_nd    DGNLperm_nd
#define SSYMperm_nd    SGNLperm_nd
#define ZSYMperm_nd    ZGNLperm_nd
#define CSYMperm_nd    CGNLperm_nd
#define ZHERperm_nd    ZGNLperm_nd
#define CHERperm_nd    CGNLperm_nd


#define DSPDperm_amf   DGNLperm_amf
#define SSPDperm_amf   SGNLperm_amf
#define ZHPDperm_amf   ZGNLperm_amf
#define CHPDperm_amf   CGNLperm_amf

#define DSPDpermamf   DGNLperm_amf
#define SSPDpermamf   SGNLperm_amf
#define ZHPDpermamf   ZGNLperm_amf
#define CHPDpermamf   CGNLperm_amf

#define DSYMperm_amf   DGNLperm_amf
#define SSYMperm_amf   SGNLperm_amf
#define ZSYMperm_amf   ZGNLperm_amf
#define CSYMperm_amf   CGNLperm_amf
#define ZHERperm_amf   ZGNLperm_amf
#define CHERperm_amf   CGNLperm_amf


#define DSPDperm_camf   DGNLperm_camf
#define SSPDperm_camf   SGNLperm_camf
#define ZHPDperm_camf   ZGNLperm_camf
#define CHPDperm_camf   CGNLperm_camf

#define DSYMperm_camf   DGNLperm_camf
#define SSYMperm_camf   SGNLperm_camf
#define ZSYMperm_camf   ZGNLperm_camf
#define CSYMperm_camf   CGNLperm_camf
#define ZHERperm_camf   ZGNLperm_camf
#define CHERperm_camf   CGNLperm_camf


#define DSPDperm_amd   DGNLperm_amd
#define SSPDperm_amd   SGNLperm_amd
#define ZHPDperm_amd   ZGNLperm_amd
#define CHPDperm_amd   CGNLperm_amd

#define DSPDperm_camd   DGNLperm_camd
#define SSPDperm_camd   SGNLperm_camd
#define ZHPDperm_camd   ZGNLperm_camd
#define CHPDperm_camd   CGNLperm_camd

#define DSPDpermamd   DGNLperm_amd
#define SSPDpermamd   SGNLperm_amd
#define ZHPDpermamd   ZGNLperm_amd
#define CHPDpermamd   CGNLperm_amd

#define DSPDpermcamd   DGNLperm_camd
#define SSPDpermcamd   SGNLperm_camd
#define ZHPDpermcamd   ZGNLperm_camd
#define CHPDpermcamd   CGNLperm_camd

#define DSYMperm_amd   DGNLperm_amd
#define SSYMperm_amd   SGNLperm_amd
#define ZSYMperm_amd   ZGNLperm_amd
#define CSYMperm_amd   CGNLperm_amd
#define ZHERperm_amd   ZGNLperm_amd
#define CHERperm_amd   CGNLperm_amd

#define DSYMperm_camd   DGNLperm_camd
#define SSYMperm_camd   SGNLperm_camd
#define ZSYMperm_camd   ZGNLperm_camd
#define CSYMperm_camd   CGNLperm_camd
#define ZHERperm_camd   ZGNLperm_camd
#define CHERperm_camd   CGNLperm_camd


#define DSPDperm_metis_e   DGNLperm_metis_e
#define SSPDperm_metis_e   SGNLperm_metis_e
#define ZHPDperm_metis_e   ZGNLperm_metis_e
#define CHPDperm_metis_e   CGNLperm_metis_e

#define DSPDperm_mtmetis   DGNLperm_mtmetis
#define SSPDperm_mtmetis   SGNLperm_mtmetis
#define ZHPDperm_mtmetis   ZGNLperm_mtmetis
#define CHPDperm_mtmetis   CGNLperm_mtmetis

#define DSPDpermmetis_e   DGNLperm_metis_e
#define SSPDpermmetis_e   SGNLperm_metis_e
#define ZHPDpermmetis_e   ZGNLperm_metis_e
#define CHPDpermmetis_e   CGNLperm_metis_e

#define DSYMperm_metis_e   DGNLperm_metis_e
#define SSYMperm_metis_e   SGNLperm_metis_e
#define ZSYMperm_metis_e   ZGNLperm_metis_e
#define CSYMperm_metis_e   CGNLperm_metis_e
#define ZHERperm_metis_e   ZGNLperm_metis_e
#define CHERperm_metis_e   CGNLperm_metis_e

#define DSYMperm_mtmetis   DGNLperm_mtmetis
#define SSYMperm_mtmetis   SGNLperm_mtmetis
#define ZSYMperm_mtmetis   ZGNLperm_mtmetis
#define CSYMperm_mtmetis   CGNLperm_mtmetis
#define ZHERperm_mtmetis   ZGNLperm_mtmetis
#define CHERperm_mtmetis   CGNLperm_mtmetis


#define DSPDperm_cmetis_e   DGNLperm_cmetis_e
#define SSPDperm_cmetis_e   SGNLperm_cmetis_e
#define ZHPDperm_cmetis_e   ZGNLperm_cmetis_e
#define CHPDperm_cmetis_e   CGNLperm_cmetis_e

#define DSYMperm_cmetis_e   DGNLperm_cmetis_e
#define SSYMperm_cmetis_e   SGNLperm_cmetis_e
#define ZSYMperm_cmetis_e   ZGNLperm_cmetis_e
#define CSYMperm_cmetis_e   CGNLperm_cmetis_e
#define ZHERperm_cmetis_e   ZGNLperm_cmetis_e
#define CHERperm_cmetis_e   CGNLperm_cmetis_e





#define DSPDperm_metis_n   DGNLperm_metis_n
#define SSPDperm_metis_n   SGNLperm_metis_n
#define ZHPDperm_metis_n   ZGNLperm_metis_n
#define CHPDperm_metis_n   CGNLperm_metis_n

#define DSPDpartition_metis_n   DGNLpartition_metis_n
#define SSPDpartition_metis_n   SGNLpartition_metis_n
#define ZHPDpartition_metis_n   ZGNLpartition_metis_n
#define CHPDpartition_metis_n   CGNLpartition_metis_n

#define DSPDpermmetis_n   DGNLperm_metis_n
#define SSPDpermmetis_n   SGNLperm_metis_n
#define ZHPDpermmetis_n   ZGNLperm_metis_n
#define CHPDpermmetis_n   CGNLperm_metis_n

#define DSYMperm_metis_n   DGNLperm_metis_n
#define SSYMperm_metis_n   SGNLperm_metis_n
#define ZSYMperm_metis_n   ZGNLperm_metis_n
#define CSYMperm_metis_n   CGNLperm_metis_n
#define ZHERperm_metis_n   ZGNLperm_metis_n
#define CHERperm_metis_n   CGNLperm_metis_n

#define DSYMpartition_metis_n   DGNLpartition_metis_n
#define SSYMpartition_metis_n   SGNLpartition_metis_n
#define ZSYMpartition_metis_n   ZGNLpartition_metis_n
#define CSYMpartition_metis_n   CGNLpartition_metis_n
#define ZHERpartition_metis_n   ZGNLpartition_metis_n
#define CHERpartition_metis_n   CGNLpartition_metis_n


#define DSPDperm_cmetis_n   DGNLperm_cmetis_n
#define SSPDperm_cmetis_n   SGNLperm_cmetis_n
#define ZHPDperm_cmetis_n   ZGNLperm_cmetis_n
#define CHPDperm_cmetis_n   CGNLperm_cmetis_n

#define DSYMperm_cmetis_n   DGNLperm_cmetis_n
#define SSYMperm_cmetis_n   SGNLperm_cmetis_n
#define ZSYMperm_cmetis_n   ZGNLperm_cmetis_n
#define CSYMperm_cmetis_n   CGNLperm_cmetis_n
#define ZHERperm_cmetis_n   ZGNLperm_cmetis_n
#define CHERperm_cmetis_n   CGNLperm_cmetis_n




#define DSPDperm_fc    DSYMperm_fc
#define SSPDperm_fc    SSYMperm_fc
#define ZHPDperm_fc    ZSYMperm_fc
#define CHPDperm_fc    CSYMperm_fc

#define DSPDperm_fcv    DSYMperm_fcv
#define SSPDperm_fcv    SSYMperm_fcv
#define ZHPDperm_fcv    ZSYMperm_fcv
#define CHPDperm_fcv    CSYMperm_fcv

#define DSPDpermfc    DSYMperm_fc
#define SSPDpermfc    SSYMperm_fc
#define ZHPDpermfc    ZSYMperm_fc
#define CHPDpermfc    CSYMperm_fc

#define DSSMperm_fc    DSYMperm_fc
#define SSSMperm_fc    SSYMperm_fc
#define ZSSMperm_fc    ZSYMperm_fc
#define CSSMperm_fc    CSYMperm_fc

#define DSSMpermfc    DSYMperm_fc
#define SSSMpermfc    SSYMperm_fc
#define ZSSMpermfc    ZSYMperm_fc
#define CSSMpermfc    CSYMperm_fc

#define ZHERperm_fc    ZSYMperm_fc
#define CHERperm_fc    CSYMperm_fc

#define ZHERpermfc    ZSYMperm_fc
#define CHERpermfc    CSYMperm_fc

#define ZSHRperm_fc    ZSYMperm_fc
#define CSHRperm_fc    CSYMperm_fc

#define ZSHRpermfc    ZSYMperm_fc
#define CSHRpermfc    CSYMperm_fc

#define DSPDperm_amd_fc    DSYMperm_amd_fc
#define SSPDperm_amd_fc    SSYMperm_amd_fc
#define ZHPDperm_amd_fc    ZSYMperm_amd_fc
#define CHPDperm_amd_fc    CSYMperm_amd_fc

#define DSPDperm_camd_fc    DSYMperm_camd_fc
#define SSPDperm_camd_fc    SSYMperm_camd_fc
#define ZHPDperm_camd_fc    ZSYMperm_camd_fc
#define CHPDperm_camd_fc    CSYMperm_camd_fc

#define DSPDperm_amd_fcv    DSYMperm_amd_fcv
#define SSPDperm_amd_fcv    SSYMperm_amd_fcv
#define ZHPDperm_amd_fcv    ZSYMperm_amd_fcv
#define CHPDperm_amd_fcv    CSYMperm_amd_fcv

#define DSPDperm_camd_fcv    DSYMperm_camd_fcv
#define SSPDperm_camd_fcv    SSYMperm_camd_fcv
#define ZHPDperm_camd_fcv    ZSYMperm_camd_fcv
#define CHPDperm_camd_fcv    CSYMperm_camd_fcv


#define DSPDperm_amf_fc    DSYMperm_amf_fc
#define SSPDperm_amf_fc    SSYMperm_amf_fc
#define ZHPDperm_amf_fc    ZSYMperm_amf_fc
#define CHPDperm_amf_fc    CSYMperm_amf_fc

#define DSPDperm_amf_fcv    DSYMperm_amf_fcv
#define SSPDperm_amf_fcv    SSYMperm_amf_fcv
#define ZHPDperm_amf_fcv    ZSYMperm_amf_fcv
#define CHPDperm_amf_fcv    CSYMperm_amf_fcv


#define DSPDperm_camf_fc    DSYMperm_camf_fc
#define SSPDperm_camf_fc    SSYMperm_camf_fc
#define ZHPDperm_camf_fc    ZSYMperm_camf_fc
#define CHPDperm_camf_fc    CSYMperm_camf_fc

#define DSPDperm_camf_fcv    DSYMperm_camf_fcv
#define SSPDperm_camf_fcv    SSYMperm_camf_fcv
#define ZHPDperm_camf_fcv    ZSYMperm_camf_fcv
#define CHPDperm_camf_fcv    CSYMperm_camf_fcv


#define DSPDperm_mmd_fc    DSYMperm_mmd_fc
#define SSPDperm_mmd_fc    SSYMperm_mmd_fc
#define ZHPDperm_mmd_fc    ZSYMperm_mmd_fc
#define CHPDperm_mmd_fc    CSYMperm_mmd_fc

#define DSPDperm_mmd_fcv    DSYMperm_mmd_fcv
#define SSPDperm_mmd_fcv    SSYMperm_mmd_fcv
#define ZHPDperm_mmd_fcv    ZSYMperm_mmd_fcv
#define CHPDperm_mmd_fcv    CSYMperm_mmd_fcv


#define DSPDperm_rcm_fc    DSYMperm_rcm_fc
#define SSPDperm_rcm_fc    SSYMperm_rcm_fc
#define ZHPDperm_rcm_fc    ZSYMperm_rcm_fc
#define CHPDperm_rcm_fc    CSYMperm_rcm_fc

#define DSPDperm_rcm_fcv    DSYMperm_rcm_fcv
#define SSPDperm_rcm_fcv    SSYMperm_rcm_fcv
#define ZHPDperm_rcm_fcv    ZSYMperm_rcm_fcv
#define CHPDperm_rcm_fcv    CSYMperm_rcm_fcv


#define DSPDperm_metis_n_fc    DSYMperm_metis_n_fc
#define SSPDperm_metis_n_fc    SSYMperm_metis_n_fc
#define ZHPDperm_metis_n_fc    ZSYMperm_metis_n_fc
#define CHPDperm_metis_n_fc    CSYMperm_metis_n_fc

#define DSPDperm_metis_n_fcv    DSYMperm_metis_n_fcv
#define SSPDperm_metis_n_fcv    SSYMperm_metis_n_fcv
#define ZHPDperm_metis_n_fcv    ZSYMperm_metis_n_fcv
#define CHPDperm_metis_n_fcv    CSYMperm_metis_n_fcv


#define DSPDperm_cmetis_n_fc    DSYMperm_cmetis_n_fc
#define SSPDperm_cmetis_n_fc    SSYMperm_cmetis_n_fc
#define ZHPDperm_cmetis_n_fc    ZSYMperm_cmetis_n_fc
#define CHPDperm_cmetis_n_fc    CSYMperm_cmetis_n_fc

#define DSPDperm_cmetis_n_fcv    DSYMperm_cmetis_n_fcv
#define SSPDperm_cmetis_n_fcv    SSYMperm_cmetis_n_fcv
#define ZHPDperm_cmetis_n_fcv    ZSYMperm_cmetis_n_fcv
#define CHPDperm_cmetis_n_fcv    CSYMperm_cmetis_n_fcv


#define DSPDperm_metis_e_fc    DSYMperm_metis_e_fc
#define SSPDperm_metis_e_fc    SSYMperm_metis_e_fc
#define ZHPDperm_metis_e_fc    ZSYMperm_metis_e_fc
#define CHPDperm_metis_e_fc    CSYMperm_metis_e_fc

#define DSPDperm_metis_e_fcv    DSYMperm_metis_e_fcv
#define SSPDperm_metis_e_fcv    SSYMperm_metis_e_fcv
#define ZHPDperm_metis_e_fcv    ZSYMperm_metis_e_fcv
#define CHPDperm_metis_e_fcv    CSYMperm_metis_e_fcv


#define DSPDperm_cmetis_e_fc    DSYMperm_cmetis_e_fc
#define SSPDperm_cmetis_e_fc    SSYMperm_cmetis_e_fc
#define ZHPDperm_cmetis_e_fc    ZSYMperm_cmetis_e_fc
#define CHPDperm_cmetis_e_fc    CSYMperm_cmetis_e_fc

#define DSPDperm_cmetis_e_fcv    DSYMperm_cmetis_e_fcv
#define SSPDperm_cmetis_e_fcv    SSYMperm_cmetis_e_fcv
#define ZHPDperm_cmetis_e_fcv    ZSYMperm_cmetis_e_fcv
#define CHPDperm_cmetis_e_fcv    CSYMperm_cmetis_e_fcv


#define DSPDperm_nd_fc    DSYMperm_nd_fc
#define SSPDperm_nd_fc    SSYMperm_nd_fc
#define ZHPDperm_nd_fc    ZSYMperm_nd_fc
#define CHPDperm_nd_fc    CSYMperm_nd_fc

#define DSPDperm_nd_fcv    DSYMperm_nd_fcv
#define SSPDperm_nd_fcv    SSYMperm_nd_fcv
#define ZHPDperm_nd_fcv    ZSYMperm_nd_fcv
#define CHPDperm_nd_fcv    CSYMperm_nd_fcv




_CPP_PREFIX integer    DSPDperm_rcm   (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
				       integer *, DILUPACKparam *);

_CPP_PREFIX integer    SSPDperm_rcm   (Smat, real *, real *, integer *,integer *,
				       integer *, SILUPACKparam *);

_CPP_PREFIX integer    ZHPDperm_rcm   (Zmat, doublecomplex *,   doublecomplex *,   integer *,integer *, 
				       integer *, ZILUPACKparam *);

_CPP_PREFIX integer    CHPDperm_rcm   (Cmat, complex *,   complex *,   integer *,integer *, 
				       integer *, CILUPACKparam *);


_CPP_PREFIX integer    DSYMperm_rcm_fc (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
					integer *, DILUPACKparam *);

_CPP_PREFIX integer    SSYMperm_rcm_fc (Smat, real *, real *, integer *,integer *,
					integer *, SILUPACKparam *);

_CPP_PREFIX integer    ZSYMperm_rcm_fc (Zmat, doublecomplex *,   doublecomplex *,   integer *,integer *, 
					integer *, ZILUPACKparam *);

_CPP_PREFIX integer    CSYMperm_rcm_fc (Cmat, complex *,   complex *,   integer *,integer *, 
					integer *, CILUPACKparam *);


_CPP_PREFIX integer    DSYMperm_rcm_fcv(Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
					integer *, DILUPACKparam *);

_CPP_PREFIX integer    SSYMperm_rcm_fcv(Smat, real *, real *, integer *,integer *,
					integer *, SILUPACKparam *);

_CPP_PREFIX integer    ZSYMperm_rcm_fcv(Zmat, doublecomplex *,   doublecomplex *,   integer *,integer *, 
					integer *, ZILUPACKparam *);

_CPP_PREFIX integer    CSYMperm_rcm_fcv(Cmat, complex *,   complex *,   integer *,integer *, 
					integer *, CILUPACKparam *);


#define DSPDpermrcm    DSPDperm_rcm   
#define SSPDpermrcm    SSPDperm_rcm   
#define ZHPDpermrcm    ZHPDperm_rcm   
#define CHPDpermrcm    CHPDperm_rcm   

#define DSYMperm_rcm   DSPDperm_rcm
#define SSYMperm_rcm   SSPDperm_rcm
#define CSYMperm_rcm   CHPDperm_rcm
#define ZSYMperm_rcm   ZHPDperm_rcm
#define CHERperm_rcm   CHPDperm_rcm
#define ZHERperm_rcm   ZHPDperm_rcm



#define DSPDperm_mmd   DGNLperm_mmd
#define SSPDperm_mmd   SGNLperm_mmd
#define ZHPDperm_mmd   ZGNLperm_mmd
#define CHPDperm_mmd   CGNLperm_mmd

#define DSPDpermmmd   DGNLperm_mmd
#define SSPDpermmmd   SGNLperm_mmd
#define ZHPDpermmmd   ZGNLperm_mmd
#define CHPDpermmmd   CGNLperm_mmd

#define DSYMperm_mmd   DGNLperm_mmd
#define SSYMperm_mmd   SGNLperm_mmd
#define ZSYMperm_mmd   ZGNLperm_mmd
#define CSYMperm_mmd   CGNLperm_mmd
#define ZHERperm_mmd   ZGNLperm_mmd
#define CHERperm_mmd   CGNLperm_mmd


#define DSPDperm_indset DGNLperm_indset
#define SSPDperm_indset SGNLperm_indset
#define ZHPDperm_indset ZGNLperm_indset
#define CHPDperm_indset CGNLperm_indset

#define DSPDpermindset DGNLperm_indset
#define SSPDpermindset SGNLperm_indset
#define ZHPDpermindset ZGNLperm_indset
#define CHPDpermindset CGNLperm_indset

#define DSYMperm_indset DGNLperm_indset
#define SSYMperm_indset SGNLperm_indset
#define ZSYMperm_indset ZGNLperm_indset
#define CSYMperm_indset CGNLperm_indset
#define ZHERperm_indset ZGNLperm_indset
#define CHERperm_indset CGNLperm_indset


_CPP_PREFIX integer    DSPDperm_pp    (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
				       integer *, DILUPACKparam *);

_CPP_PREFIX integer    DSPDperm_cpp   (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
				       integer *, DILUPACKparam *);

_CPP_PREFIX integer    SSPDperm_pp    (Smat, real *, real *, integer *,integer *,
				       integer *, SILUPACKparam *);
_CPP_PREFIX integer    SSPDperm_cpp   (Smat, real *, real *, integer *,integer *,
				       integer *, SILUPACKparam *);

_CPP_PREFIX integer    ZHPDperm_pp    (Zmat, doublecomplex *,   doublecomplex *, integer *,integer *,
				       integer *, ZILUPACKparam *);

_CPP_PREFIX integer    ZHPDperm_cpp   (Zmat, doublecomplex *,   doublecomplex *, integer *,integer *,
				       integer *, ZILUPACKparam *);

_CPP_PREFIX integer    CHPDperm_pp    (Cmat, complex *,   complex *, integer *,integer *,
				       integer *, CILUPACKparam *);

_CPP_PREFIX integer    CHPDperm_cpp    (Cmat, complex *,   complex *, integer *,integer *,
				       integer *, CILUPACKparam *);

#define DSPDpermpp     DSPDperm_pp    
#define SSPDpermpp     SSPDperm_pp    
#define ZHPDpermpp     ZHPDperm_pp    
#define CHPDpermpp     CHPDperm_pp    
#define DSPDpermcpp     DSPDperm_cpp    
#define SSPDpermcpp     SSPDperm_cpp    
#define ZHPDpermcpp     ZHPDperm_cpp    
#define CHPDpermcpp     CHPDperm_cpp    




#define DGNLpermnull      DGNLperm_null   
#define DGNLpermcnull      DGNLperm_cnull   
#define DGNLpermnd        DGNLperm_nd     
#define DGNLpermrcm       DGNLperm_rcm    
#define DGNLpermamf       DGNLperm_amf    
#define DGNLpermcamf       DGNLperm_camf    
#define DGNLpermamd       DGNLperm_amd 
#define DGNLpermcamd      DGNLperm_camd 
#define DGNLpermmmd       DGNLperm_mmd    
#define DGNLpermpq        DGNLperm_pq     
#define DGNLpermfc        DGNLperm_fc     
#define DSYMpermfc        DSYMperm_fc     
#define DGNLpermp         DGNLperm_p      
#define DGNLpermindset    DGNLperm_indset         


#define DGNLpermmwm_rcm          DGNLperm_mwm_rcm       
#define DGNLpermmwm_mmd          DGNLperm_mwm_mmd       
#define DGNLpermmwm_amf          DGNLperm_mwm_amf       
#define DGNLpermmwm_camf          DGNLperm_mwm_camf       
#define DGNLpermmwm_amd          DGNLperm_mwm_amd
#define DGNLpermmwm_camd         DGNLperm_mwm_camd
#define DGNLpermmwm_metis_e      DGNLperm_mwm_metis_e    
#define DGNLpermmwm_cmetis_e      DGNLperm_mwm_cmetis_e    
#define DGNLpermmwm_metis_n      DGNLperm_mwm_metis_n    
#define DGNLpermmwm_cmetis_n      DGNLperm_mwm_cmetis_n    

#define DGNLpermmatching_rcm          DGNLperm_matching_rcm       
#define DGNLpermmatching_mmd          DGNLperm_matching_mmd       
#define DGNLpermmatching_amf          DGNLperm_matching_amf       
#define DGNLpermmatching_camf          DGNLperm_matching_camf       
#define DGNLpermmatching_amd          DGNLperm_matching_amd
#define DGNLpermmatching_camd         DGNLperm_matching_camd
#define DGNLpermmatching_metis_e      DGNLperm_matching_metis_e    
#define DGNLpermmatching_cmetis_e      DGNLperm_matching_cmetis_e    
#define DGNLpermmatching_metis_n      DGNLperm_matching_metis_n    
#define DGNLpermmatching_cmetis_n      DGNLperm_matching_cmetis_n    
                                                       
#define DGNLpermmc64_rcm         DGNLperm_mc64_rcm      
#define DGNLpermmc64_mtmetis     DGNLperm_mc64_mtmetis
#define DGNLpermmc64_mmd         DGNLperm_mc64_mmd      
#define DGNLpermmc64_amf         DGNLperm_mc64_amf      
#define DGNLpermmc64_camf         DGNLperm_mc64_camf      
#define DGNLpermmc64_amd         DGNLperm_mc64_amd
#define DGNLpermmc64_camd        DGNLperm_mc64_camd
#define DGNLpermmc64_metis_e     DGNLperm_mc64_metis_e  
#define DGNLpermmc64_cmetis_e     DGNLperm_mc64_cmetis_e  
#define DGNLpermmc64_metis_n     DGNLperm_mc64_metis_n  
#define DGNLpermmc64_cmetis_n     DGNLperm_mc64_cmetis_n  


#define SGNLpermnull    SGNLperm_null    
#define SGNLpermcnull    SGNLperm_cnull    
#define SGNLpermnd      SGNLperm_nd      
#define SGNLpermrcm     SGNLperm_rcm     
#define SGNLpermamf     SGNLperm_amf     
#define SGNLpermcamf     SGNLperm_camf     
#define SGNLpermamd     SGNLperm_amd
#define SGNLpermcamd    SGNLperm_camd
#define SGNLpermmmd     SGNLperm_mmd     
#define SGNLpermpq      SGNLperm_pq      
#define SGNLpermfc      SGNLperm_fc      
#define SSYMpermfc      SSYMperm_fc      
#define SGNLpermp       SGNLperm_p       
#define SGNLpermindset  SGNLperm_indset  

#define SSPDpermrcm	    SSPDperm_rcm	
#define SSPDpermpp	    SSPDperm_pp	
#define SSPDpermcpp	    SSPDperm_cpp	

#define SGNLpermmwm_rcm       SGNLperm_mwm_rcm    
#define SGNLpermmwm_mmd       SGNLperm_mwm_mmd    
#define SGNLpermmwm_amf       SGNLperm_mwm_amf    
#define SGNLpermmwm_camf       SGNLperm_mwm_camf    
#define SGNLpermmwm_amd       SGNLperm_mwm_amd
#define SGNLpermmwm_camd      SGNLperm_mwm_camd
#define SGNLpermmwm_metis_e   SGNLperm_mwm_metis_e
#define SGNLpermmwm_cmetis_e   SGNLperm_mwm_cmetis_e
#define SGNLpermmwm_metis_n   SGNLperm_mwm_metis_n
#define SGNLpermmwm_cmetis_n   SGNLperm_mwm_cmetis_n

#define SGNLpermmatching_rcm       SGNLperm_matching_rcm    
#define SGNLpermmatching_mmd       SGNLperm_matching_mmd    
#define SGNLpermmatching_amf       SGNLperm_matching_amf    
#define SGNLpermmatching_camf       SGNLperm_matching_camf    
#define SGNLpermmatching_amd       SGNLperm_matching_amd
#define SGNLpermmatching_camd      SGNLperm_matching_camd
#define SGNLpermmatching_metis_e   SGNLperm_matching_metis_e
#define SGNLpermmatching_cmetis_e   SGNLperm_matching_cmetis_e
#define SGNLpermmatching_metis_n   SGNLperm_matching_metis_n
#define SGNLpermmatching_cmetis_n   SGNLperm_matching_cmetis_n


#define SGNLpermmc64_rcm          SGNLperm_mc64_rcm    
#define SGNLpermmc64_mtmetis      SGNLperm_mc64_metis
#define SGNLpermmc64_mmd	  SGNLperm_mc64_mmd    
#define SGNLpermmc64_amf	  SGNLperm_mc64_amf    
#define SGNLpermmc64_camf	  SGNLperm_mc64_camf    
#define SGNLpermmc64_amd	  SGNLperm_mc64_amd
#define SGNLpermmc64_camd	  SGNLperm_mc64_camd
#define SGNLpermmc64_metis_e	  SGNLperm_mc64_metis_e
#define SGNLpermmc64_cmetis_e	  SGNLperm_mc64_cmetis_e
#define SGNLpermmc64_metis_n	  SGNLperm_mc64_metis_n
#define SGNLpermmc64_cmetis_n	  SGNLperm_mc64_cmetis_n


#define CGNLpermnull       CGNLperm_null  
#define CGNLpermcnull       CGNLperm_cnull  
#define CGNLpermnd	   CGNLperm_nd    
#define CGNLpermrcm	   CGNLperm_rcm   
#define CGNLpermamf	   CGNLperm_amf   
#define CGNLpermcamf	   CGNLperm_camf   
#define CGNLpermamd	   CGNLperm_amd
#define CGNLpermcamd	   CGNLperm_camd
#define CGNLpermmmd	   CGNLperm_mmd   
#define CGNLpermpq	   CGNLperm_pq    
#define CGNLpermfc	   CGNLperm_fc    
#define CSYMpermfc	   CSYMperm_fc    
#define CGNLpermp	   CGNLperm_p     
#define CGNLpermindset	   CGNLperm_indset

#define CGNLpermmwm_rcm      CGNLperm_mwm_rcm	   
#define CGNLpermmwm_mmd	     CGNLperm_mwm_mmd	   
#define CGNLpermmwm_amf	     CGNLperm_mwm_amf	   
#define CGNLpermmwm_camf     CGNLperm_mwm_camf	   
#define CGNLpermmwm_amd	     CGNLperm_mwm_amd
#define CGNLpermmwm_camd     CGNLperm_mwm_camd
#define CGNLpermmwm_metis_e  CGNLperm_mwm_metis_e
#define CGNLpermmwm_cmetis_e  CGNLperm_mwm_cmetis_e
#define CGNLpermmwm_metis_n  CGNLperm_mwm_metis_n
#define CGNLpermmwm_cmetis_n  CGNLperm_mwm_cmetis_n

#define CGNLpermmatching_rcm      CGNLperm_matching_rcm	   
#define CGNLpermmatching_mmd	     CGNLperm_matching_mmd	   
#define CGNLpermmatching_amf	     CGNLperm_matching_amf	   
#define CGNLpermmatching_camf	     CGNLperm_matching_camf	   
#define CGNLpermmatching_amd	     CGNLperm_matching_amd
#define CGNLpermmatching_camd	     CGNLperm_matching_camd
#define CGNLpermmatching_metis_e  CGNLperm_matching_metis_e
#define CGNLpermmatching_cmetis_e  CGNLperm_matching_cmetis_e
#define CGNLpermmatching_metis_n  CGNLperm_matching_metis_n
#define CGNLpermmatching_cmetis_n  CGNLperm_matching_cmetis_n


#define CGNLpermmc64_rcm      CGNLperm_mc64_rcm    
#define CGNLpermmc64_mtmetis  CGNLperm_mc64_mtmetis
#define CGNLpermmc64_mmd      CGNLperm_mc64_mmd    
#define CGNLpermmc64_amf      CGNLperm_mc64_amf    
#define CGNLpermmc64_camf      CGNLperm_mc64_camf    
#define CGNLpermmc64_amd      CGNLperm_mc64_amd 
#define CGNLpermmc64_camd     CGNLperm_mc64_camd 
#define CGNLpermmc64_metis_e  CGNLperm_mc64_metis_e
#define CGNLpermmc64_cmetis_e  CGNLperm_mc64_cmetis_e
#define CGNLpermmc64_metis_n  CGNLperm_mc64_metis_n
#define CGNLpermmc64_cmetis_n  CGNLperm_mc64_cmetis_n

#define CHPDpermrcm	 CHPDperm_rcm   
#define CHPDpermpp	 CHPDperm_pp    
#define CHPDpermcpp	 CHPDperm_cpp    


#define ZGNLpermnull      ZGNLperm_null	 
#define ZGNLpermcnull      ZGNLperm_cnull	 
#define ZGNLpermnd	  ZGNLperm_nd	 
#define ZGNLpermrcm	  ZGNLperm_rcm	 
#define ZGNLpermamf	  ZGNLperm_amf	 
#define ZGNLpermcamf	  ZGNLperm_camf	 
#define ZGNLpermamd	  ZGNLperm_amd
#define ZGNLpermcamd	  ZGNLperm_camd
#define ZGNLpermmmd	  ZGNLperm_mmd	 
#define ZGNLpermpq	  ZGNLperm_pq	 
#define ZGNLpermfc	  ZGNLperm_fc	 
#define ZSYMpermfc	  ZSYMperm_fc	 
#define ZGNLpermp	  ZGNLperm_p	 
#define ZGNLpermindset    ZGNLperm_indset   

#define ZGNLpermmwm_rcm         ZGNLperm_mwm_rcm	       
#define ZGNLpermmwm_mmd		ZGNLperm_mwm_mmd	       
#define ZGNLpermmwm_amf		ZGNLperm_mwm_amf	       
#define ZGNLpermmwm_camf	ZGNLperm_mwm_camf	       
#define ZGNLpermmwm_amd		ZGNLperm_mwm_amd
#define ZGNLpermmwm_camd	ZGNLperm_mwm_camd
#define ZGNLpermmwm_metis_e	ZGNLperm_mwm_metis_e    
#define ZGNLpermmwm_cmetis_e	ZGNLperm_mwm_cmetis_e    
#define ZGNLpermmwm_metis_n	ZGNLperm_mwm_metis_n    
#define ZGNLpermmwm_cmetis_n	ZGNLperm_mwm_cmetis_n    
				                       
#define ZGNLpermmatching_rcm         ZGNLperm_matching_rcm	       
#define ZGNLpermmatching_mmd		ZGNLperm_matching_mmd	       
#define ZGNLpermmatching_amf		ZGNLperm_matching_amf	       
#define ZGNLpermmatching_camf		ZGNLperm_matching_camf	       
#define ZGNLpermmatching_amd		ZGNLperm_matching_amd
#define ZGNLpermmatching_camd		ZGNLperm_matching_camd
#define ZGNLpermmatching_metis_e	ZGNLperm_matching_metis_e    
#define ZGNLpermmatching_cmetis_e	ZGNLperm_matching_cmetis_e    
#define ZGNLpermmatching_metis_n	ZGNLperm_matching_metis_n    
#define ZGNLpermmatching_cmetis_n	ZGNLperm_matching_cmetis_n    
				                       
#define ZGNLpermmc64_rcm	ZGNLperm_mc64_rcm       
#define ZGNLpermmc64_mtmetis	ZGNLperm_mc64_mtmetis       
#define ZGNLpermmc64_mmd	ZGNLperm_mc64_mmd       
#define ZGNLpermmc64_amf	ZGNLperm_mc64_amf       
#define ZGNLpermmc64_camf	ZGNLperm_mc64_camf       
#define ZGNLpermmc64_amd	ZGNLperm_mc64_amd
#define ZGNLpermmc64_camd	ZGNLperm_mc64_camd
#define ZGNLpermmc64_metis_e	ZGNLperm_mc64_metis_e   
#define ZGNLpermmc64_cmetis_e	ZGNLperm_mc64_cmetis_e   
#define ZGNLpermmc64_metis_n	ZGNLperm_mc64_metis_n   
#define ZGNLpermmc64_cmetis_n	ZGNLperm_mc64_cmetis_n   
				                       
#define ZHPDpermrcm		ZHPDperm_rcm	       
#define ZHPDpermpp		ZHPDperm_pp	       
#define ZHPDpermcpp		ZHPDperm_cpp	       


#define ZHERindfc  ZSYMindfc               
#define CHERindfc  CSYMindfc               



_CPP_PREFIX void swapj(integer *, integer, integer);
_CPP_PREFIX void dswapm(doubleprecision *, integer, integer);
_CPP_PREFIX void sswapm(real *, integer, integer);
_CPP_PREFIX integer DPQpermF(Dmat, integer, integer *, integer *, integer *, doubleprecision,doubleprecision *,    integer *);
_CPP_PREFIX integer ZPQpermF(Zmat, integer, integer *, integer *, integer *, doubleprecision,doublecomplex *, integer *);
_CPP_PREFIX integer indAMF(Dmat, integer, integer *, integer *, doubleprecision);

_CPP_PREFIX integer Dindset(Dmat, integer, integer *, integer *, doubleprecision);
_CPP_PREFIX integer Sindset(Smat, integer, integer *, integer *, real);
_CPP_PREFIX integer Zindset(Zmat, integer, integer *, integer *, doubleprecision);
_CPP_PREFIX integer Cindset(Cmat, integer, integer *, integer *, real);

_CPP_PREFIX void Dindfc(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void Sindfc(Smat, integer *, integer *, integer *, real, real *, integer *);
_CPP_PREFIX void Zindfc(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void Cindfc(Cmat, integer *, integer *, integer *, real, real *, integer *);

_CPP_PREFIX void Dindfc_rs(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void Sindfc_rs(Smat, integer *, integer *, integer *, real, real *, integer *);
_CPP_PREFIX void Zindfc_rs(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void Cindfc_rs(Cmat, integer *, integer *, integer *, real, real *, integer *);

_CPP_PREFIX void DSYMindfc(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void SSYMindfc(Smat, integer *, integer *, integer *, real, real *, integer *);
_CPP_PREFIX void ZSYMindfc(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void CSYMindfc(Cmat, integer *, integer *, integer *, real, real *, integer *);

_CPP_PREFIX void DSYMindfc_rs(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void SSYMindfc_rs(Smat, integer *, integer *, integer *, real, real *, integer *);
_CPP_PREFIX void ZSYMindfc_rs(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX void CSYMindfc_rs(Cmat, integer *, integer *, integer *, real, real *, integer *);


_CPP_PREFIX void SSYMpindfc(Smat, integer *,integer *, integer *,integer *,
			    integer *, real, real *,integer *);
_CPP_PREFIX void DSYMpindfc(Dmat, integer *,integer *, integer *,integer *,
			    integer *, doubleprecision, doubleprecision *,integer *);
_CPP_PREFIX void CSYMpindfc(Cmat, integer *,integer *, integer *,integer *,
			    integer *, real, real *,integer *);
_CPP_PREFIX void ZSYMpindfc(Zmat, integer *,integer *, integer *,integer *,
			    integer *, doubleprecision, doubleprecision *,integer *);
_CPP_PREFIX void CHERpindfc(Cmat, integer *,integer *, integer *,integer *,
			    integer *, real, real *,integer *);
_CPP_PREFIX void ZHERpindfc(Zmat, integer *,integer *, integer *,integer *,
			    integer *, doubleprecision, doubleprecision *,integer *);

_CPP_PREFIX void SSYMpindfc_rs(Smat, integer *,integer *, integer *,integer *,
			       integer *, real, real *,integer *);
_CPP_PREFIX void DSYMpindfc_rs(Dmat, integer *,integer *, integer *,integer *,
			       integer *, doubleprecision, doubleprecision *,integer *);
_CPP_PREFIX void CSYMpindfc_rs(Cmat, integer *,integer *, integer *,integer *,
			       integer *, real, real *,integer *);
_CPP_PREFIX void ZSYMpindfc_rs(Zmat, integer *,integer *, integer *,integer *,
			       integer *, doubleprecision, doubleprecision *,integer *);
_CPP_PREFIX void CHERpindfc_rs(Cmat, integer *,integer *, integer *,integer *,
			       integer *, real, real *,integer *);
_CPP_PREFIX void ZHERpindfc_rs(Zmat, integer *,integer *, integer *,integer *,
			       integer *, doubleprecision, doubleprecision *,integer *);


_CPP_PREFIX void SSYMpindfcv(Smat, integer *,integer *, integer *,integer *,
			     integer *, real, real *, integer, 
			     real *,integer *,real *);
_CPP_PREFIX void DSYMpindfcv(Dmat, integer *,integer *, integer *,integer *,
			     integer *, doubleprecision, doubleprecision *, integer, 
			     doubleprecision *,integer *,doubleprecision *);
_CPP_PREFIX void CSYMpindfcv(Cmat, integer *,integer *, integer *,integer *,
			     integer *, real, complex *, integer, 
			     real *,integer *,complex *);
_CPP_PREFIX void ZSYMpindfcv(Zmat, integer *,integer *, integer *,integer *,
			     integer *, doubleprecision, doublecomplex *, integer, 
			     doubleprecision *,integer *,doublecomplex *);
_CPP_PREFIX void CHERpindfcv(Cmat, integer *,integer *, integer *,integer *,
			     integer *, real, complex *, integer,
			     real *,integer *,complex *);
_CPP_PREFIX void ZHERpindfcv(Zmat, integer *,integer *, integer *,integer *,
			     integer *, doubleprecision, doublecomplex *, integer, 
			     doubleprecision *,integer *,doublecomplex *);

_CPP_PREFIX void SSYMpindfcv_rs(Smat, integer *,integer *, integer *,integer *,
				integer *, real, real *, integer, 
				real *,integer *,real *);
_CPP_PREFIX void DSYMpindfcv_rs(Dmat, integer *,integer *, integer *,integer *,
				integer *, doubleprecision, doubleprecision *, integer, 
				doubleprecision *,integer *,doubleprecision *);
_CPP_PREFIX void CSYMpindfcv_rs(Cmat, integer *,integer *, integer *,integer *,
				integer *, real, complex *, integer, 
				real *,integer *,complex *);
_CPP_PREFIX void ZSYMpindfcv_rs(Zmat, integer *,integer *, integer *,integer *,
				integer *, doubleprecision, doublecomplex *, integer, 
				doubleprecision *,integer *,doublecomplex *);
_CPP_PREFIX void CHERpindfcv_rs(Cmat, integer *,integer *, integer *,integer *,
				integer *, real, complex *, integer,
				real *,integer *,complex *);
_CPP_PREFIX void ZHERpindfcv_rs(Zmat, integer *,integer *, integer *,integer *,
				integer *, doubleprecision, doublecomplex *, integer, 
				doubleprecision *,integer *,doublecomplex *);


_CPP_PREFIX void SSYMbuildblock(Smat, integer *,integer *, integer *, real *);
_CPP_PREFIX void DSYMbuildblock(Dmat, integer *,integer *, integer *, doubleprecision *);
_CPP_PREFIX void CSYMbuildblock(Cmat, integer *,integer *, integer *, real *);
_CPP_PREFIX void ZSYMbuildblock(Zmat, integer *,integer *, integer *, doubleprecision *);
_CPP_PREFIX void CHERbuildblock(Cmat, integer *,integer *, integer *, real *);
_CPP_PREFIX void ZHERbuildblock(Zmat, integer *,integer *, integer *, doubleprecision *);







_CPP_PREFIX void Dindfcv(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer, doubleprecision *, integer *);
_CPP_PREFIX void Sindfcv(Smat, integer *, integer *, integer *, real,            real *,            integer, real *,            integer *);
_CPP_PREFIX void Zindfcv(Zmat, integer *, integer *, integer *, doubleprecision, doublecomplex *,   integer, doubleprecision *, integer *);
_CPP_PREFIX void Cindfcv(Cmat, integer *, integer *, integer *, real,            complex *,         integer, real *,            integer *);

_CPP_PREFIX void DSYMindfcv(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer, doubleprecision *, integer *);
_CPP_PREFIX void SSYMindfcv(Smat, integer *, integer *, integer *, real,            real *,            integer, real *,            integer *);
_CPP_PREFIX void ZSYMindfcv(Zmat, integer *, integer *, integer *, doubleprecision, doublecomplex *,   integer, doubleprecision *, integer *);
_CPP_PREFIX void CSYMindfcv(Cmat, integer *, integer *, integer *, real,            complex *,         integer, real *,            integer *);
_CPP_PREFIX void ZHERindfcv(Zmat, integer *, integer *, integer *, doubleprecision, doublecomplex *,   integer, doubleprecision *, integer *);
_CPP_PREFIX void CHERindfcv(Cmat, integer *, integer *, integer *, real,            complex *,         integer, real *,            integer *);


_CPP_PREFIX void dqsortr2i(doubleprecision *, integer *, integer *, integer, integer);
_CPP_PREFIX void sqsortr2i(real *, integer *, integer *, integer, integer);

_CPP_PREFIX void Dclear(integer, doubleprecision *, integer);
_CPP_PREFIX void Sclear(integer, real *,            integer);
_CPP_PREFIX void Zclear(integer, doublecomplex *,   integer);
_CPP_PREFIX void Cclear(integer, complex *,         integer);

_CPP_PREFIX void IP_etree(integer *, integer *, integer, integer *, 
			  integer *, integer *, integer *);
_CPP_PREFIX void Setree(SSparseMatrix *, integer *, 
			integer *, integer *, integer *);
_CPP_PREFIX void Detree(DSparseMatrix *, integer *, 
			integer *, integer *, integer *);
_CPP_PREFIX void Cetree(CSparseMatrix *, integer *, 
			integer *, integer *, integer *);
_CPP_PREFIX void Zetree(ZSparseMatrix *, integer *, 
			integer *, integer *, integer *);
_CPP_PREFIX void Setree0(SSparseMatrix *, integer *, 
			 integer *);
_CPP_PREFIX void Detree0(DSparseMatrix *, integer *, 
			 integer *);
_CPP_PREFIX void Cetree0(CSparseMatrix *, integer *, 
			 integer *);
_CPP_PREFIX void Zetree0(ZSparseMatrix *, integer *, 
			 integer *);
_CPP_PREFIX void IP_post_order(integer *, integer, integer *, integer *);
_CPP_PREFIX integer IP_tdfs(integer, integer, integer *, integer *, 
			    integer *, integer *);


/* ********************************************* */
/* ******      Definitions for solvers     ***** */
_CPP_PREFIX void      Dpcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      Dfpcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      Dfpcgnosave(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *, 
				    integer *,doubleprecision *,doubleprecision *,integer *);
_CPP_PREFIX void      Dbcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      DSYMbcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);

_CPP_PREFIX void      DSYMqmr(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      DSYMfqmr(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      DSYMfqmrnosave(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *,
				       integer *,doubleprecision *,doubleprecision *,integer *);

_CPP_PREFIX void      Dgmres(integer *,doubleprecision *,doubleprecision *,integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      Dfgmres(integer *,doubleprecision *,doubleprecision *,integer *,doubleprecision *,doubleprecision *);
_CPP_PREFIX void      Dfgmresnosave(integer *,doubleprecision *,doubleprecision *,integer *,doubleprecision *,doubleprecision *,
				      integer *,doubleprecision *,doubleprecision *,integer *);

_CPP_PREFIX void      Spcg(integer *,real *,   real *,   integer *,real *,real *);
_CPP_PREFIX void      Sfpcg(integer *,real *,   real *,   integer *,real *,real *);
_CPP_PREFIX void      Sfpcgnosave(integer *,real *,   real *,   integer *,real *,real *,
				    integer *,real *,real *,integer *);
_CPP_PREFIX void      Sbcg(integer *,real *,   real *,   integer *,real *,real *);
_CPP_PREFIX void      SSYMbcg(integer *,real *,   real *,   integer *,real *,real *);

_CPP_PREFIX void      SSYMqmr(integer *,real *,   real *,   integer *,real *,real *);
_CPP_PREFIX void      SSYMfqmr(integer *,real *,   real *,   integer *,real *,real *);
_CPP_PREFIX void      SSYMfqmrnosave(integer *,real *,   real *,   integer *,real *,real *,
				       integer *,real *,real *,integer *);

_CPP_PREFIX void      Sgmres(integer *,real *,real *,integer *,real *,real *);
_CPP_PREFIX void      Sfgmres(integer *,real *,real *,integer *,real *,real *);
_CPP_PREFIX void      Sfgmresnosave(integer *,real *,real *,integer *,real *,real *,
				      integer *,real *,real *,integer *);



_CPP_PREFIX void      Zpcg(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      Zfpcg(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      Zfpcgnosave(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *,
				    integer *,doubleprecision *,doublecomplex *,integer *);
_CPP_PREFIX void      Zbcg(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZSYMbcg(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZHERbcg(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);

_CPP_PREFIX void      ZSYMqmr(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZSYMfqmr(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZSYMfqmrnosave(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *,
				       integer *,doubleprecision *,doublecomplex *,integer *);
_CPP_PREFIX void      ZHERqmr(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZHERfqmr(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      ZHERfqmrnosave(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *,
				       integer *,doubleprecision *,doublecomplex *,integer *);

_CPP_PREFIX void      Zgmres(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      Zfgmres(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *);
_CPP_PREFIX void      Zfgmresnosave(integer *,doublecomplex *,doublecomplex *,integer *,doubleprecision *,doublecomplex *,
				      integer *,doubleprecision *,doublecomplex *,integer *);

_CPP_PREFIX void      Cpcg(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      Cfpcg(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      Cfpcgnosave(integer *,complex *,complex *,integer *,real *,complex *,
				    integer *,real *,complex *,integer *);
_CPP_PREFIX void      Cbcg(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CSYMbcg(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CHERbcg(integer *,complex *,complex *,integer *,real *,complex *);

_CPP_PREFIX void      CSYMqmr(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CSYMfqmr(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CSYMfqmrnosave(integer *,complex *,complex *,integer *,real *,complex *,
				       integer *,real *,complex *,integer *);
_CPP_PREFIX void      CHERqmr(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CHERfqmr(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      CHERfqmrnosave(integer *,complex *,complex *,integer *,real *,complex *,
				       integer *,real *,complex *,integer *);			       

_CPP_PREFIX void      Cgmres(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      Cfgmres(integer *,complex *,complex *,integer *,real *,complex *);
_CPP_PREFIX void      Cfgmresnosave(integer *,complex *,complex *,integer *,real *,complex *,
				      integer *,real *,complex *,integer *);


_CPP_PREFIX doubleprecision Ddistdot(integer *,doubleprecision *,integer *,doubleprecision *,integer *);

_CPP_PREFIX real            Sdistdot(integer *,real *,integer *,real *,integer *);

_CPP_PREFIX doublecomplex   Zdistdotc(integer *,doublecomplex *,integer *,doublecomplex *,integer *);
_CPP_PREFIX doublecomplex   Zdistdotu(integer *,doublecomplex *,integer *,doublecomplex *,integer *);

_CPP_PREFIX complex         Cdistdotc(integer *,complex *,integer *,complex *,integer *);
_CPP_PREFIX complex         Cdistdotu(integer *,complex *,integer *,complex *,integer *);

_CPP_PREFIX integer AMGsolver(SPARSEmat *, AMGlevelmat *, ILUPACKparam *, 
			      void *, void *);

_CPP_PREFIX integer DGNLAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				  doubleprecision *, doubleprecision *);
_CPP_PREFIX integer DGNLSYMAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				     doubleprecision *, doubleprecision *);
_CPP_PREFIX integer DGNLSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				     doubleprecision *, doubleprecision *);
_CPP_PREFIX integer DSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				  doubleprecision *, doubleprecision *);
_CPP_PREFIX integer DSYMAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				  doubleprecision *, doubleprecision *);
_CPP_PREFIX integer DSYMSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
				     doubleprecision *, doubleprecision *);


_CPP_PREFIX integer SGNLAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				  real *, real *);
_CPP_PREFIX integer SGNLSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				     real *, real *);
_CPP_PREFIX integer SGNLSYMAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				     real *, real *);
_CPP_PREFIX integer SSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				  real *, real *);
_CPP_PREFIX integer SSYMAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				  real *, real *);
_CPP_PREFIX integer SSYMSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
				     real *, real *);

_CPP_PREFIX integer ZGNLAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				  doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLDGNLAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				     doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLHERAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				     doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLSYMAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				     doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZGNLDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				  doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHPDDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHERHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				     doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHERDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHERDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZHERAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				  doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZSYMAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
				  doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZSYMDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);
_CPP_PREFIX integer ZSYMDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
				      doublecomplex *, doublecomplex *);

_CPP_PREFIX integer CGNLAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				  complex *, complex *);
_CPP_PREFIX integer CGNLSGNLAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CGNLSYMAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				     complex *, complex *);
_CPP_PREFIX integer CGNLHERAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				     complex *, complex *);
_CPP_PREFIX integer CGNLHPDAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				     complex *, complex *);
_CPP_PREFIX integer CGNLSSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CGNLSSYMAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CHPDAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				  complex *, complex *);
_CPP_PREFIX integer CHPDSSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CHERAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				  complex *, complex *);
_CPP_PREFIX integer CHERHPDAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				     complex *, complex *);
_CPP_PREFIX integer CHERSSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CHERSSYMAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CSYMAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
				  complex *, complex *);
_CPP_PREFIX integer CSYMSSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);
_CPP_PREFIX integer CSYMSSYMAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
				      complex *, complex *);

_CPP_PREFIX void AMGinit(SPARSEmat *, ILUPACKparam *);

_CPP_PREFIX void DGNLAMGinit(Dmat *, DILUPACKparam *);
_CPP_PREFIX void DGNLAMGgetparams(DILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *, integer *);
_CPP_PREFIX void DGNLAMGsetparams(Dmat *, DILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer , integer );

_CPP_PREFIX void SGNLAMGinit(Smat *, SILUPACKparam *);
_CPP_PREFIX void SGNLAMGgetparams(SILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *, integer *);
_CPP_PREFIX void SGNLAMGsetparams(Smat *, SILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer , integer );

_CPP_PREFIX void ZGNLAMGinit(Zmat *, ZILUPACKparam *);
_CPP_PREFIX void ZGNLAMGgetparams(ZILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *, integer *);
_CPP_PREFIX void ZGNLAMGsetparams(Zmat *, ZILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer , integer );

_CPP_PREFIX void CGNLAMGinit(Cmat *, CILUPACKparam *);
_CPP_PREFIX void CGNLAMGgetparams(CILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *, integer *);
_CPP_PREFIX void CGNLAMGsetparams(Cmat *, CILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer , integer );


_CPP_PREFIX void DSPDAMGinit(Dmat *, DILUPACKparam *);
_CPP_PREFIX void DSPDAMGgetparams(DILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *);
_CPP_PREFIX void DSPDAMGsetparams(Dmat *, DILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer );

_CPP_PREFIX void DSYMAMGinit(Dmat *, DILUPACKparam *);
_CPP_PREFIX void DSYMAMGgetparams(DILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *);
_CPP_PREFIX void DSYMAMGsetparams(Dmat *, DILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer );

_CPP_PREFIX void SSPDAMGinit(Smat *, SILUPACKparam *);
_CPP_PREFIX void SSPDAMGgetparams(SILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *);
_CPP_PREFIX void SSPDAMGsetparams(Smat *, SILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer );

_CPP_PREFIX void SSYMAMGinit(Smat *, SILUPACKparam *);
_CPP_PREFIX void SSYMAMGgetparams(SILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *);
_CPP_PREFIX void SSYMAMGsetparams(Smat *, SILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer );

_CPP_PREFIX void ZHPDAMGinit(Zmat *, ZILUPACKparam *);
_CPP_PREFIX void ZHERAMGinit(Zmat *, ZILUPACKparam *);
_CPP_PREFIX void ZSYMAMGinit(Zmat *, ZILUPACKparam *);
_CPP_PREFIX void ZHPDAMGgetparams(ZILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *);
_CPP_PREFIX void ZHPDAMGsetparams(Zmat *, ZILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer );
_CPP_PREFIX void ZHERAMGgetparams(ZILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *);
_CPP_PREFIX void ZHERAMGsetparams(Zmat *, ZILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer );
_CPP_PREFIX void ZSYMAMGgetparams(ZILUPACKparam *,
				  integer *, integer *, doubleprecision *, doubleprecision *,
				  doubleprecision *, integer *);
_CPP_PREFIX void ZSYMAMGsetparams(Zmat *, ZILUPACKparam *,
				  integer , integer , doubleprecision *, doubleprecision ,
				  doubleprecision , integer );

_CPP_PREFIX void CHPDAMGinit(Cmat *, CILUPACKparam *);
_CPP_PREFIX void CHERAMGinit(Cmat *, CILUPACKparam *);
_CPP_PREFIX void CSYMAMGinit(Cmat *, CILUPACKparam *);
_CPP_PREFIX void CHPDAMGgetparams(CILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *);
_CPP_PREFIX void CHPDAMGsetparams(Cmat *, CILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer );
_CPP_PREFIX void CHERAMGgetparams(CILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *);
_CPP_PREFIX void CHERAMGsetparams(Cmat *, CILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer );
_CPP_PREFIX void CSYMAMGgetparams(CILUPACKparam *,
				  integer *, integer *, real *, real *,
				  real *, integer *);
_CPP_PREFIX void CSYMAMGsetparams(Cmat *, CILUPACKparam *,
				  integer , integer , real *, real ,
				  real , integer );
_CPP_PREFIX void Sconsistent_threads(Smat, SILUPACKparam *);
_CPP_PREFIX void Dconsistent_threads(Dmat, DILUPACKparam *);
_CPP_PREFIX void Cconsistent_threads(Cmat, CILUPACKparam *);
_CPP_PREFIX void Zconsistent_threads(Zmat, ZILUPACKparam *);



/* ********************************************* */
/* ******      Definitions for sparskit    ***** */
_CPP_PREFIX void      Dcsrcsc(integer *, integer *, integer *, doubleprecision *, integer *, integer *, 
			      doubleprecision *, integer *, integer *);
_CPP_PREFIX void      Dcsrcsc(integer *,integer *,integer *,doubleprecision *,
			      integer *,integer *,doubleprecision *,integer *,integer *);

_CPP_PREFIX void      Scsrcsc(integer *, integer *, integer *, real *, integer *, integer *, 
			      real *, integer *, integer *);
_CPP_PREFIX void      Scsrcsc(integer *,integer *,integer *,real *,
			      integer *,integer *,real *,integer *,integer *);

_CPP_PREFIX void      Zcsrcsc(integer *, integer *, integer *, doublecomplex *, integer *, integer *, 
			      doublecomplex *, integer *, integer *);
_CPP_PREFIX void      Zcsrcsc (integer *,integer *,integer *,doublecomplex *,
			       integer *,integer *,doublecomplex *,integer *,integer *);

_CPP_PREFIX void      Ccsrcsc(integer *, integer *, integer *, complex *, integer *, integer *, 
			      complex *, integer *, integer *);
_CPP_PREFIX void      Ccsrcsc(integer *,integer *,integer *,complex *,
			      integer *,integer *,complex *,integer *,integer *);




/* ******************************************* */
/* ******      Definitions for tools     ***** */
_CPP_PREFIX void SPARSEmatvec(SPARSEmat, void *, void *);
_CPP_PREFIX void DGNLmatvec(Dmat, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLmattvec(Dmat, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DGNLmathvec(Dmat, doubleprecision *, doubleprecision *);

_CPP_PREFIX void SGNLmatvec(Smat, real *, real *);
_CPP_PREFIX void SGNLmattvec(Smat, real *, real *);
_CPP_PREFIX void SGNLmathvec(Smat, real *, real *);

_CPP_PREFIX void ZGNLmatvec(Zmat, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZGNLmattvec(Zmat, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZGNLmathvec(Zmat, doublecomplex *, doublecomplex *);

_CPP_PREFIX void CGNLmatvec(Cmat, complex *, complex *);
_CPP_PREFIX void CGNLmattvec(Cmat, complex *, complex *);
_CPP_PREFIX void CGNLmathvec(Cmat, complex *, complex *);


_CPP_PREFIX void DSYMmatvec(Dmat, doubleprecision *, doubleprecision *);
_CPP_PREFIX void DSSMmatvec(Dmat, doubleprecision *, doubleprecision *);

_CPP_PREFIX void SSYMmatvec(Smat, real *, real *);
_CPP_PREFIX void SSSMmatvec(Smat, real *, real *);

_CPP_PREFIX void ZHERmatvec(Zmat, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZSHRmatvec(Zmat, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZSYMmatvec(Zmat, doublecomplex *, doublecomplex *);
_CPP_PREFIX void ZSSMmatvec(Zmat, doublecomplex *, doublecomplex *);

_CPP_PREFIX void CHERmatvec(Cmat, complex *, complex *);
_CPP_PREFIX void CSHRmatvec(Cmat, complex *, complex *);
_CPP_PREFIX void CSYMmatvec(Cmat, complex *, complex *);
_CPP_PREFIX void CSSMmatvec(Cmat, complex *, complex *);


_CPP_PREFIX void Sqsort(real *,            integer *, integer *, integer *); 
_CPP_PREFIX void Dqsort(doubleprecision *, integer *, integer *, integer *); 
_CPP_PREFIX void Cqsort(complex *,         integer *, integer *, integer *); 
_CPP_PREFIX void Zqsort(doublecomplex *,   integer *, integer *, integer *); 
_CPP_PREFIX void Sqsortpair(real *,            real *,           integer *, integer *, integer *); 
_CPP_PREFIX void Dqsortpair(doubleprecision *, doubleprecision *,integer *, integer *, integer *); 
_CPP_PREFIX void Cqsortpair(complex *,         complex *,        integer *, integer *, integer *); 
_CPP_PREFIX void Zqsortpair(doublecomplex *,   doublecomplex *,  integer *, integer *, integer *); 
_CPP_PREFIX void Sqsort2(real *,            integer *, integer *, integer *); 
_CPP_PREFIX void Dqsort2(doubleprecision *, integer *, integer *, integer *); 
_CPP_PREFIX void Cqsort2(complex *,         integer *, integer *, integer *); 
_CPP_PREFIX void Zqsort2(doublecomplex *,   integer *, integer *, integer *); 

_CPP_PREFIX void Sbqsort(real *,           integer *,integer *,integer *,integer *); 
_CPP_PREFIX void Dbqsort(doubleprecision *,integer *,integer *,integer *,integer *); 
_CPP_PREFIX void Cbqsort(complex *,        integer *,integer *,integer *,integer *); 
_CPP_PREFIX void Zbqsort(doublecomplex *,  integer *,integer *,integer *,integer *); 


_CPP_PREFIX void qqsorti(integer *, integer *, integer *); 

_CPP_PREFIX integer  Dspartran(Dmat, Dmat *, integer, integer);
_CPP_PREFIX void Dsetupgraph(Dmat, Dmat *, integer *, integer *, size_t);
_CPP_PREFIX void Dsetupcgraph(Dmat, Dmat *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Dsetupgraph_epsilon(Dmat, Dmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t);
_CPP_PREFIX void Dsetupgraph_epsilon_sp(Dmat, Dmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Dqqsort(doubleprecision *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Dqqsort2(doubleprecision *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Dqqsortgnl(doubleprecision *, integer *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Dqsortgnl(doubleprecision *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Dqqsorts(doubleprecision *, integer *, integer *, integer *); 
_CPP_PREFIX void Dqqsorts2(doubleprecision *, integer *, integer *, integer *); 
_CPP_PREFIX void Dcperm(Dmat *, integer *);
_CPP_PREFIX void Drperm(Dmat *, integer *);

_CPP_PREFIX integer  Sspartran(Smat, Smat *, integer, integer);
_CPP_PREFIX void Ssetupgraph(Smat, Smat *, integer *, integer *, size_t);
_CPP_PREFIX void Ssetupcgraph(Smat, Smat *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Ssetupgraph_epsilon(Smat, Smat *, real, real *, integer *, integer *, size_t);
_CPP_PREFIX void Ssetupgraph_epsilon_sp(Smat, Smat *, real, real *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Sqqsort(real *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Sqqsortgnl(real *, integer *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Sqsortgnl(real *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Scperm(Smat *, integer *);
_CPP_PREFIX void Srperm(Smat *, integer *);

_CPP_PREFIX integer  Zspartran(Zmat, Zmat *, integer, integer);
_CPP_PREFIX void Zsetupgraph(Zmat, Zmat *, integer *, integer *, size_t);
_CPP_PREFIX void Zsetupcgraph(Zmat, Zmat *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Zsetupgraph_epsilon(Zmat, Zmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t);
_CPP_PREFIX void Zsetupgraph_epsilon_sp(Zmat, Zmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Zqqsort(doublecomplex *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Zqqsortgnl(doublecomplex *, integer *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Zqsortgnl(doublecomplex *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Zcperm(Zmat *, integer *);
_CPP_PREFIX void Zrperm(Zmat *, integer *);

_CPP_PREFIX integer  Cspartran(Cmat, Cmat *, integer, integer);
_CPP_PREFIX void Csetupgraph(Cmat, Cmat *, integer *, integer *, size_t);
_CPP_PREFIX void Csetupcgraph(Cmat, Cmat *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Csetupgraph_epsilon(Cmat, Cmat *, real, real *, integer *, integer *, size_t);
_CPP_PREFIX void Csetupgraph_epsilon_sp(Cmat, Cmat *, real, real *, integer *, integer *, size_t, integer *);
_CPP_PREFIX void Cqqsort(complex *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Cqqsortgnl(complex *, integer *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Cqsortgnl(complex *, integer *, integer *, integer *, integer *); 
_CPP_PREFIX void Ccperm(Cmat *, integer *);
_CPP_PREFIX void Crperm(Cmat *, integer *);


_CPP_PREFIX void *MAlloc(size_t, char *);
_CPP_PREFIX void *CAlloc(size_t, size_t, char *);
_CPP_PREFIX void *ReAlloc(void *, size_t, char *);
_CPP_PREFIX void *FRee(void *);
_CPP_PREFIX doubleprecision dgeteps();
_CPP_PREFIX real            sgeteps();
_CPP_PREFIX float evaluate_time(float *, float *);
_CPP_PREFIX float evaluatetime(float *, float *);


_CPP_PREFIX void Droscal(integer *,integer *,integer *,doubleprecision *,integer *,
			 integer *,doubleprecision *,doubleprecision *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Dcoscal(integer *,integer *,integer *,doubleprecision *,integer *,
			 integer *,doubleprecision *,doubleprecision *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Drowscale(integer *,integer *,doubleprecision *,integer *,
			   integer *,doubleprecision *,integer *);
_CPP_PREFIX void Dcolscale(integer *,integer *,doubleprecision *,integer *,
			   integer *,doubleprecision *,integer *);

_CPP_PREFIX void Drowcscale(integer *,integer *,doubleprecision *,integer *,
			    integer *,doubleprecision *,
			    integer *,doubleprecision *,integer *);
_CPP_PREFIX void Dcolcscale(integer *,integer *,doubleprecision *,integer *,
			    integer *,doubleprecision *,
			    integer *,doubleprecision *,integer *);
_CPP_PREFIX void Drowcpscale(integer *,integer *,doubleprecision *,integer *,
			     integer *,doubleprecision *,
			     integer *,integer *);
_CPP_PREFIX void Dcolcpscale(integer *,integer *,doubleprecision *,integer *,
			     integer *,doubleprecision *,
			     integer *,integer *);
_CPP_PREFIX void DSYMcscale(integer *, doubleprecision *,integer *,
			    integer *,doubleprecision *,doubleprecision *,
			    integer *,doubleprecision *,integer *);
_CPP_PREFIX void DSYMcpscale(integer *, doubleprecision *,integer *,
			     integer *,doubleprecision *,doubleprecision *,
			     integer *,integer *);

_CPP_PREFIX void DSPDscale  (integer *, doubleprecision *,integer *,
			     integer *,doubleprecision *,integer *);
_CPP_PREFIX void DSPDcscale (integer *, doubleprecision *,integer *,
			     integer *,doubleprecision *,integer *,doubleprecision *,integer *);
_CPP_PREFIX void DSPDcpscale(integer *, doubleprecision *,integer *,
			     integer *,doubleprecision *,integer *,integer *);

_CPP_PREFIX void DSYMscale(integer *, doubleprecision *,integer *,
			   integer *,doubleprecision *,doubleprecision *,
			   integer *);

_CPP_PREFIX void Sroscal(integer *,integer *,integer *,real *,integer *,
			 integer *,real *,real *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Scoscal(integer *,integer *,integer *,real *,integer *,
			 integer *,real *,real *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Srowscale(integer *,integer *,real *,integer *,
			   integer *,real *,integer *);
_CPP_PREFIX void Scolscale(integer *,integer *,real *,integer *,
			   integer *,real *,integer *);

_CPP_PREFIX void Srowcscale(integer *,integer *,real *,integer *,
			    integer *,real *,
			    integer *,real *,integer *);
_CPP_PREFIX void Scolcscale(integer *,integer *,real *,integer *,
			    integer *,real *,
			    integer *,real *,integer *);
_CPP_PREFIX void Srowcpscale(integer *,integer *,real *,integer *,
			     integer *,real *,
			     integer *,integer *);
_CPP_PREFIX void Scolcpscale(integer *,integer *,real *,integer *,
			     integer *,real *,
			     integer *,integer *);
_CPP_PREFIX void SSYMcscale(integer *, real *,integer *,
			    integer *,real *,real *,
			    integer *,real *,integer *);
_CPP_PREFIX void SSYMcpscale(integer *, real *,integer *,
			     integer *,real *,real *,
			     integer *,integer *);

_CPP_PREFIX void SSPDscale(integer *, real *,integer *,
			   integer *,real *,integer *);
_CPP_PREFIX void SSPDcscale(integer *, real *,integer *,
			   integer *,real *,integer *,real *,integer *);
_CPP_PREFIX void SSPDcpscale(integer *, real *,integer *,
			   integer *,real *,integer *,integer *);

_CPP_PREFIX void SSYMscale(integer *, real *,integer *,
			   integer *,real *,real *,integer *);

_CPP_PREFIX void Zroscal(integer *,integer *,integer *,doublecomplex *,integer *,
			 integer *,doublecomplex *,doublecomplex *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Zcoscal(integer *,integer *,integer *,doublecomplex *,integer *,
			 integer *,doublecomplex *,doublecomplex *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Zrowscale(integer *,integer *,doublecomplex *,integer *,
			   integer *,doublecomplex *,integer *);
_CPP_PREFIX void Zcolscale(integer *,integer *,doublecomplex *,integer *,
			   integer *,doublecomplex *,integer *);

_CPP_PREFIX void Zrowcscale(integer *,integer *,doublecomplex *,integer *,
			    integer *,doublecomplex *,
			    integer *,doublecomplex *,integer *);
_CPP_PREFIX void Zcolcscale(integer *,integer *,doublecomplex *,integer *,
			    integer *,doublecomplex *,
			    integer *,doublecomplex *,integer *);
_CPP_PREFIX void Zrowcpscale(integer *,integer *,doublecomplex *,integer *,
			     integer *,doublecomplex *,
			     integer *,integer *);
_CPP_PREFIX void Zcolcpscale(integer *,integer *,doublecomplex *,integer *,
			     integer *,doublecomplex *,
			     integer *,integer *);
_CPP_PREFIX void ZSYMcscale(integer *, doublecomplex *,integer *,
			    integer *,doublecomplex *,doublecomplex *,
			    integer *,doublecomplex *,integer *);
_CPP_PREFIX void ZHERcscale(integer *, doublecomplex *,integer *,
			    integer *,doublecomplex *,doublecomplex *,
			    integer *,doublecomplex *,integer *);
_CPP_PREFIX void ZSYMcpscale(integer *, doublecomplex *,integer *,
			     integer *,doublecomplex *,doublecomplex *,
			     integer *,integer *);
_CPP_PREFIX void ZHERcpscale(integer *, doublecomplex *,integer *,
			     integer *,doublecomplex *,doublecomplex *,
			     integer *,integer *);

_CPP_PREFIX void ZHPDscale  (integer *, doublecomplex *,integer *,
			     integer *,doublecomplex *,integer *);
_CPP_PREFIX void ZHPDcscale (integer *, doublecomplex *,integer *,
			     integer *,doublecomplex *,integer *,doublecomplex *,integer *);
_CPP_PREFIX void ZHPDcpscale(integer *, doublecomplex *,integer *,
			     integer *,doublecomplex *,integer *,integer *);

_CPP_PREFIX void ZSYMscale(integer *, doublecomplex *,integer *,
			   integer *,doublecomplex *,doublecomplex *,integer *);
_CPP_PREFIX void ZHERscale(integer *, doublecomplex *,integer *,
			   integer *,doublecomplex *,doublecomplex *,integer *);

_CPP_PREFIX void Croscal(integer *,integer *,integer *,complex *,integer *,
			 integer *,complex *,complex *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Ccoscal(integer *,integer *,integer *,complex *,integer *,
			 integer *,complex *,complex *,integer *,integer *,
			 integer *);
_CPP_PREFIX void Crowscale(integer *,integer *,complex *,integer *,
			   integer *,complex *,integer *);
_CPP_PREFIX void Ccolscale(integer *,integer *,complex *,integer *,
			   integer *,complex *,integer *);

_CPP_PREFIX void Crowcscale(integer *,integer *,complex *,integer *,
			    integer *,complex *,
			    integer *,complex *,integer *);
_CPP_PREFIX void Ccolcscale(integer *,integer *,complex *,integer *,
			    integer *,complex *,
			    integer *,complex *,integer *);
_CPP_PREFIX void Crowcpscale(integer *,integer *,complex *,integer *,
			     integer *,complex *,
			     integer *,integer *);
_CPP_PREFIX void Ccolcpscale(integer *,integer *,complex *,integer *,
			     integer *,complex *,
			     integer *,integer *);
_CPP_PREFIX void CSYMcscale(integer *, complex *,integer *,
			    integer *, complex *,complex *,
			    integer *, complex *,integer *);
_CPP_PREFIX void CHERcscale(integer *, complex *,integer *,
			    integer *,complex *,complex *,
			    integer *, complex *,integer *);
_CPP_PREFIX void CSYMcpscale(integer *, complex *,integer *,
			     integer *, complex *,complex *,
			     integer *, integer *);
_CPP_PREFIX void CHERcpscale(integer *, complex *,integer *,
			     integer *,complex *,complex *,
			     integer *, integer *);

_CPP_PREFIX void CHPDscale  (integer *, complex *,integer *,
			     integer *,complex *,integer *);
_CPP_PREFIX void CHPDcscale (integer *, complex *,integer *,
			     integer *,complex *,integer *,complex *,integer *);
_CPP_PREFIX void CHPDcpscale (integer *, complex *,integer *,
			      integer *,complex *,integer *,integer *);

_CPP_PREFIX void CSYMscale(integer *, complex *,integer *,
			   integer *, complex *,complex *,integer *);
_CPP_PREFIX void CHERscale(integer *, complex *,integer *,
			   integer *,complex *,complex *,integer *);


_CPP_PREFIX integer DPQpermF(Dmat, integer, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX integer SPQpermF(Smat, integer, integer *, integer *, integer *, real,            real *,            integer *);
_CPP_PREFIX integer ZPQpermF(Zmat, integer, integer *, integer *, integer *, doubleprecision, doublecomplex *,   integer *);
_CPP_PREFIX integer CPQpermF(Cmat, integer, integer *, integer *, integer *, real,            complex *,         integer *);

_CPP_PREFIX integer DPPpermF(Dmat, integer, integer *, integer *, doubleprecision, doubleprecision *, integer *);
_CPP_PREFIX integer SPPpermF(Smat, integer, integer *, integer *, real,            real *,            integer *);
_CPP_PREFIX integer ZPPpermF(Zmat, integer, integer *, integer *, doubleprecision, doublecomplex *,   integer *);
_CPP_PREFIX integer CPPpermF(Cmat, integer, integer *, integer *, real,            complex *,         integer *);





_CPP_PREFIX void      Dreadmtc(integer *,integer *,integer *,character *,doubleprecision *,integer *,
			       integer *,doubleprecision *,integer *,character *,integer *,integer *,
			       integer *,character *,character *,character *,
			       integer *, integer *, integer *, doubleprecision *,
			       integer *,ftnlen,ftnlen, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Sreadmtc(integer *,integer *,integer *,character *,real *,integer *,
			       integer *,real *,integer *,character *,integer *,integer *,
			       integer *,character *,character *,character *,
			       integer *, integer *, integer *, real *,
			       integer *,ftnlen,ftnlen, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Zreadmtc(integer *,integer *,integer *,character *,doublecomplex *,integer *,
			       integer *,doublecomplex *,integer *,character *,integer *,integer *,
			       integer *,character *,character *,character *,
			       integer *, integer *, integer *, doublecomplex *,
			       integer *,ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Creadmtc(integer *,integer *,integer *,character *,complex *,integer *,
			       integer *,complex *,integer *,character *,integer *,integer *,
			       integer *,character *,character *,character *,
			       integer *, integer *, integer *, complex *,
			       integer *,ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Dwritemtc(character *, doubleprecision *,integer *,integer *,
				doubleprecision *,integer *, character *,
				integer *,integer *,integer *,
				character *,character *,character *,
				ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Swritemtc(character *, real *,integer *,integer *,
				real *,integer *, character *,
				integer *,integer *,integer *,
				character *,character *,character *,
				ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Zwritemtc(character *, doublecomplex *,integer *,integer *,
				doublecomplex *,integer *, character *,
				integer *,integer *,integer *,
				character *,character *,character *,
				ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Cwritemtc(character *, complex *,integer *,integer *,
				complex *,integer *, character *,
				integer *,integer *,integer *,
				character *,character *,character *,
				ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);

_CPP_PREFIX void      Dreadvectors(character *,doubleprecision *,integer *,integer *,
				   character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Sreadvectors(character *,real *,integer *,integer *,
				   character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Zreadvectors(character *,doublecomplex *,integer *,integer *,
				   character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Creadvectors(character *,complex *,integer *,integer *,
				   character *,character *, ftnlen,ftnlen,ftnlen);

_CPP_PREFIX void      Dwritevectors(character *,doubleprecision *,integer *,integer *,
				    character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Swritevectors(character *,real *,integer *,integer *,
				    character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Zwritevectors(character *,doublecomplex *,integer *,integer *,
				    character *,character *, ftnlen,ftnlen,ftnlen);
_CPP_PREFIX void      Cwritevectors(character *,complex *,integer *,integer *,
				    character *,character *, ftnlen,ftnlen,ftnlen);


_CPP_PREFIX integer  DSSMsmwm(Dmat, integer *, integer *, doubleprecision *),
     SSSMsmwm(Smat, integer *, integer *, real *),
     DSYMsmwm(Dmat, integer *, integer *, doubleprecision *),
     SSYMsmwm(Smat, integer *, integer *, real *),
     CSSMsmwm(Cmat, integer *, integer *, complex *),
     ZSSMsmwm(Zmat, integer *, integer *, doublecomplex *),
     CSYMsmwm(Cmat, integer *, integer *, complex *),
     ZSYMsmwm(Zmat, integer *, integer *, doublecomplex *),
     CSHRsmwm(Cmat, integer *, integer *, complex *),
     ZSHRsmwm(Zmat, integer *, integer *, doublecomplex *),
     CHERsmwm(Cmat, integer *, integer *, complex *),
     ZHERsmwm(Zmat, integer *, integer *, doublecomplex *),
     DGNLsmwm(Dmat, integer *, integer *, doubleprecision *),
     SGNLsmwm(Smat, integer *, integer *, real *),
     CGNLsmwm(Cmat, integer *, integer *, complex *),
     ZGNLsmwm(Zmat, integer *, integer *, doublecomplex *);




_CPP_PREFIX void SSYMAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX void DSYMAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX void CSYMAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZSYMAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX void CHERAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZHERAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);


_CPP_PREFIX void SSPDAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX void DSPDAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX void CHPDAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZHPDAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
		        
		        
_CPP_PREFIX void SGNLAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX void DGNLAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX void CGNLAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZGNLAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX void AMGdelete(SPARSEmat *, AMGlevelmat *, ILUPACKparam *);

_CPP_PREFIX integer dsymilupack   (integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
				   doubleprecision *, doubleprecision *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer dsymilupackfac(long *, long *, 
				   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer dsymilupacksol(long *, long *, 
				   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer dsymilupackdel(long *, long *, 
				   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);


_CPP_PREFIX integer ssymilupack   (integer *, integer *, integer *, real *, real *, real *,
				   real *, real *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer ssymilupackfac(long *, long *, 
				   integer *, integer *, integer *, real *, real *, real *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer ssymilupacksol(long *, long *, 
				   integer *, integer *, integer *, real *, real *, real *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer ssymilupackdel(long *, long *, 
				   integer *, integer *, integer *, real *, real *, real *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);


_CPP_PREFIX integer csymilupack   (integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer csymilupackfac(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer csymilupacksol(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer csymilupackdel(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer chermilupack   (integer *, integer *, integer *, complex *, complex *, complex *,
				    real *, real *, integer *, integer *, 
				    real *, integer *);
_CPP_PREFIX integer cherilupackfac(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer cherilupacksol(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);
_CPP_PREFIX integer cherilupackdel(long *, long *, 
				   integer *, integer *, integer *, complex *, complex *, complex *,
				   real *, real *, integer *, integer *, integer *, 
				   real *, integer *);


_CPP_PREFIX integer zsymilupack   (integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zsymilupackfac(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zsymilupacksol(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zsymilupackdel(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);

_CPP_PREFIX integer zherilupack   (integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zherilupackfac(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zherilupacksol(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);
_CPP_PREFIX integer zherilupackdel(long *, long *, 
				   integer *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
				   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
				   doubleprecision *, integer *);




_CPP_PREFIX void sgnlamginit(integer *, integer *, integer *,
			     real *, integer *, character *, real *,real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void dgnlamginit(integer *, integer *, integer *,
			     double *, integer *, character *, double *,double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void cgnlamginit(integer *, integer *, integer *,
			     complex *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void zgnlamginit(integer *, integer *, integer *,
			     doublecomplex *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX void sspdamginit(integer *, integer *, integer *,
			     real *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void dspdamginit(integer *, integer *, integer *,
			     double *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void chpdamginit(integer *, integer *, integer *,
			     complex *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void zhpdamginit(integer *, integer *, integer *,
			     doublecomplex *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX void ssymamginit(integer *, integer *, integer *,
			     real *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void dsymamginit(integer *, integer *, integer *,
			     double *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void cheramginit(integer *, integer *, integer *,
			     complex *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void zheramginit(integer *, integer *, integer *,
			     doublecomplex *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void csymamginit(integer *, integer *, integer *,
			     complex *, integer *, character *, real *, real *,
			     real *, real *, integer *, real *,
			     integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX void zsymamginit(integer *, integer *, integer *,
			     doublecomplex *, integer *, character *, double *, double *,
			     double *, double *, integer *, double *,
			     integer *, integer *, integer *, integer *, integer *);


_CPP_PREFIX int sgnlamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dgnlamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int cgnlamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zgnlamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX int sspdamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dspdamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int chpdamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zhpdamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX int ssymamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dsymamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int cheramgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zheramgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int csymamgfactor(size_t *, size_t *, 
			      integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zsymamgfactor(size_t *, size_t *,
			      integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);


_CPP_PREFIX int sgnlamgsolver(size_t *, size_t *, 
			      real *, real *, integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dgnlamgsolver(size_t *, size_t *, 
			      double *, double *, integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int cgnlamgsolver(size_t *, size_t *, 
			      complex *, complex *, integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX int zgnlamgsolver(size_t *, size_t *, 
			      doublecomplex *, doublecomplex *, integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX int sspdamgsolver(size_t *, size_t *, 
			      real *, real *, integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dspdamgsolver(size_t *, size_t *, 
			      double *, double *, integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int chpdamgsolver(size_t *, size_t *, 
			      complex *, complex *, integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zhpdamgsolver(size_t *, size_t *, 
			      doublecomplex *, doublecomplex *, integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);

_CPP_PREFIX int ssymamgsolver(size_t *, size_t *, 
			      real *, real *, integer *, integer *, integer *,
			      real *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int dsymamgsolver(size_t *, size_t *, 
			      double *, double *, integer *, integer *, integer *,
			      double *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int cheramgsolver(size_t *, size_t *, 
			      complex *, complex *, integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zheramgsolver(size_t *, size_t *, 
			      doublecomplex *, doublecomplex *, integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int csymamgsolver(size_t *, size_t *, 
			      complex *, complex *, integer *, integer *, integer *,
			      complex *, integer *, character *, real *, real *,
			      real *, real *, integer *, real *,
			      integer *, integer *, integer *, integer *, integer *);
_CPP_PREFIX int zsymamgsolver(size_t *, size_t *, 
			      doublecomplex *, doublecomplex *, integer *, integer *, integer *,
			      doublecomplex *, integer *, character *, double *, double *,
			      double *, double *, integer *, double *,
			      integer *, integer *, integer *, integer *, integer *);


_CPP_PREFIX void sgnlamgsol(size_t *, size_t *, 
			    real *, real *, integer *);
_CPP_PREFIX void dgnlamgsol(size_t *, size_t *, 
			    double *, double *, integer *);
_CPP_PREFIX void cgnlamgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zgnlamgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);

_CPP_PREFIX void sspdamgsol(size_t *, size_t *, 
			    real *, real *, integer *);
_CPP_PREFIX void dspdamgsol(size_t *, size_t *, 
			    double *, double *, integer *);
_CPP_PREFIX void chpdamgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zhpdamgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);

_CPP_PREFIX void ssymamgsol(size_t *, size_t *, 
			    real *, real *, integer *);
_CPP_PREFIX void dsymamgsol(size_t *, size_t *, 
			    double *, double *, integer *);
_CPP_PREFIX void csymamgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zsymamgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);
_CPP_PREFIX void cheramgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zheramgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);
_CPP_PREFIX void ssymamgsol(size_t *, size_t *, 
			    real *, real *, integer *);
_CPP_PREFIX void dsymamgsol(size_t *, size_t *, 
			    double *, double *, integer *);
_CPP_PREFIX void csymamgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zsymamgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);
_CPP_PREFIX void cheramgsol(size_t *, size_t *, 
			    complex *, complex *, integer *);
_CPP_PREFIX void zheramgsol(size_t *, size_t *, 
			    doublecomplex *, doublecomplex *, integer *);


_CPP_PREFIX void sgnlamgdelete(size_t *, size_t *);
_CPP_PREFIX void dgnlamgdelete(size_t *, size_t *);
_CPP_PREFIX void cgnlamgdelete(size_t *, size_t *);
_CPP_PREFIX void zgnlamgdelete(size_t *, size_t *);

_CPP_PREFIX void sspdamgdelete(size_t *, size_t *);
_CPP_PREFIX void dspdamgdelete(size_t *, size_t *);
_CPP_PREFIX void chpdamgdelete(size_t *, size_t *);
_CPP_PREFIX void zhpdamgdelete(size_t *, size_t *);

_CPP_PREFIX void ssymamgdelete(size_t *, size_t *);
_CPP_PREFIX void dsymamgdelete(size_t *, size_t *);
_CPP_PREFIX void cheramgdelete(size_t *, size_t *);
_CPP_PREFIX void zheramgdelete(size_t *, size_t *);
_CPP_PREFIX void csymamgdelete(size_t *, size_t *);
_CPP_PREFIX void zsymamgdelete(size_t *, size_t *);



_CPP_PREFIX void sspdamginfo(size_t *, size_t *, integer *, integer *, 
			     integer *, real *);
_CPP_PREFIX void dspdamginfo(size_t *, size_t *, 
			     integer *, integer *, integer *, double *);
_CPP_PREFIX void chpdamginfo(size_t *, size_t *, 
			     integer *, integer *, integer *, complex *);
_CPP_PREFIX void zhpdamginfo(size_t *, size_t *, 
			     integer *, integer *, integer *, doublecomplex *);

_CPP_PREFIX void ssymamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, real *);
_CPP_PREFIX void dsymamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, double *);
_CPP_PREFIX void cheramginfo(size_t *, size_t *,
			     integer *, integer *, integer *, complex *);
_CPP_PREFIX void zheramginfo(size_t *, size_t *,
			     integer *, integer *, integer *, doublecomplex *);
_CPP_PREFIX void csymamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, complex *);
_CPP_PREFIX void zsymamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, doublecomplex *);

_CPP_PREFIX void sgnlamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, real *);
_CPP_PREFIX void dgnlamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, double *);
_CPP_PREFIX void cgnlamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, complex *);
_CPP_PREFIX void zgnlamginfo(size_t *, size_t *,
			     integer *, integer *, integer *, doublecomplex *);


_CPP_PREFIX size_t sspdamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t dspdamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t chpdamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t zhpdamgnnz(size_t *, size_t *);

_CPP_PREFIX size_t ssymamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t dsymamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t cheramgnnz(size_t *, size_t *);
_CPP_PREFIX size_t zheramgnnz(size_t *, size_t *);
_CPP_PREFIX size_t csymamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t zsymamgnnz(size_t *, size_t *);

_CPP_PREFIX size_t sgnlamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t dgnlamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t cgnlamgnnz(size_t *, size_t *);
_CPP_PREFIX size_t zgnlamgnnz(size_t *, size_t *);

_CPP_PREFIX void ssymspdamgconvert(size_t *, size_t *);
_CPP_PREFIX void dsymspdamgconvert(size_t *, size_t *);
_CPP_PREFIX void cherhpdamgconvert(size_t *, size_t *);
_CPP_PREFIX void zherhpdamgconvert(size_t *, size_t *);

_CPP_PREFIX void samgundoscaling(size_t *, size_t *, integer *, integer *, integer *, real *);
_CPP_PREFIX void damgundoscaling(size_t *, size_t *, integer *, integer *, integer *, doubleprecision *);
_CPP_PREFIX void camgundoscaling(size_t *, size_t *, integer *, integer *, integer *, complex *);
_CPP_PREFIX void zamgundoscaling(size_t *, size_t *, integer *, integer *, integer *, doublecomplex *);


/*
void DGNLSYM(Dmat , Dmat *, integer *);
void SGNLSYM(Smat , Smat *, integer *);
void CGNLSYM(Cmat , Cmat *, integer *);
void ZGNLSYM(Zmat , Zmat *, integer *);
void CGNLHER(Cmat , Cmat *, integer *);
void ZGNLHER(Zmat , Zmat *, integer *);
*/
_CPP_PREFIX void SSYMSPDAMGconvert(SAMGlevelmat *);
_CPP_PREFIX void DSYMSPDAMGconvert(DAMGlevelmat *);
_CPP_PREFIX void CHERHPDAMGconvert(CAMGlevelmat *);
_CPP_PREFIX void ZHERHPDAMGconvert(ZAMGlevelmat *);



_CPP_PREFIX void SSYMpiluclsol(integer *,real *,real *,real *, integer *);
_CPP_PREFIX void SSYMpilucusol(integer *,real *,real *,real *, integer *);
_CPP_PREFIX void DSYMpiluclsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			       integer *);
_CPP_PREFIX void DSYMpilucusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
			       integer *);
_CPP_PREFIX void CSYMpiluclsol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMpilucusol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERpiluclsol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERpilucusol(integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void ZSYMpiluclsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZSYMpilucusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZHERpiluclsol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZHERpilucusol(integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);


_CPP_PREFIX void SSYMppiluclsol(integer *,integer *,real *,real *,real *, integer *);
_CPP_PREFIX void SSYMppilucusol(integer *,integer *,real *,real *,real *, integer *);
_CPP_PREFIX void SSYMbppiluclsol(integer *,integer *,real *,real *,real *, integer *);
_CPP_PREFIX void SSYMbppilucusol(integer *,integer *,real *,real *,real *, integer *);
_CPP_PREFIX void DSYMppiluclsol(integer *,integer *,doubleprecision *,doubleprecision *,
				doubleprecision *, integer *);
_CPP_PREFIX void DSYMppilucusol(integer *,integer *,doubleprecision *,doubleprecision *,
				doubleprecision *, integer *);
_CPP_PREFIX void DSYMbppiluclsol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *, integer *);
_CPP_PREFIX void DSYMbppilucusol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *, integer *);
_CPP_PREFIX void CSYMppiluclsol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMppilucusol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMbppiluclsol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CSYMbppilucusol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERppiluclsol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERppilucusol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERbppiluclsol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void CHERbppilucusol(integer *,integer *,complex *,complex *,complex *, integer *);
_CPP_PREFIX void ZSYMppiluclsol(integer *,integer *,doublecomplex *,doublecomplex *,
				doublecomplex *, integer *);
_CPP_PREFIX void ZSYMppilucusol(integer *,integer *,doublecomplex *,doublecomplex *,
				doublecomplex *, integer *);
_CPP_PREFIX void ZSYMbppiluclsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZSYMbppilucusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZHERppiluclsol(integer *,integer *,doublecomplex *,doublecomplex *,
				doublecomplex *, integer *);
_CPP_PREFIX void ZHERppilucusol(integer *,integer *,doublecomplex *,doublecomplex *,
				doublecomplex *, integer *);
_CPP_PREFIX void ZHERbppiluclsol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);
_CPP_PREFIX void ZHERbppilucusol(integer *,integer *,doublecomplex *,doublecomplex *,doublecomplex *, integer *);


_CPP_PREFIX integer Smps_arms(Smat, integer *, real *, real *,  
			      integer *, real *);
_CPP_PREFIX integer Dmps_arms(Dmat, integer *, doubleprecision *, doubleprecision *, 
			      integer *, doubleprecision *);
_CPP_PREFIX integer Cmps_arms(Cmat, integer *, real *, real *,  
			      integer *, real *);
_CPP_PREFIX integer Zmps_arms(Zmat, integer *, doubleprecision *, doubleprecision *, 
			      integer *, doubleprecision *);

_CPP_PREFIX void sprivatesptrs(character *uplo, integer *n, integer *nrhs, real *ap, integer *ipiv, real *b, integer *ldb, integer *info, ftnlen uplolen);
_CPP_PREFIX void dprivatesptrs(character *uplo, integer *n, integer *nrhs, doubleprecision *ap, integer *ipiv, doubleprecision *b, integer *ldb, integer *info, ftnlen uplolen);
_CPP_PREFIX void cprivatehptrs(character *uplo, integer *n, integer *nrhs, complex *ap, integer *ipiv, complex *b, integer *ldb, integer *info, ftnlen uplolen);
_CPP_PREFIX void zprivatehptrs(character *uplo, integer *n, integer *nrhs, doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, integer *info, ftnlen uplolen);


_CPP_PREFIX void SGNLAMGsetupparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DGNLAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CGNLAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZGNLAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);

_CPP_PREFIX void SGNLAMGdeleteparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DGNLAMGdeleteparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CGNLAMGdeleteparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZGNLAMGdeleteparameters(Zmat *, ZILUPACKparam *, integer);

_CPP_PREFIX void SSPDAMGsetupparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DSPDAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CHPDAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZHPDAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);

_CPP_PREFIX void SSPDAMGdeleteparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DSPDAMGdeleteparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CHPDAMGdeleteparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZHPDAMGdeleteparameters(Zmat *, ZILUPACKparam *, integer);

_CPP_PREFIX void SSYMAMGsetupparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DSYMAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CHERAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZHERAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);
_CPP_PREFIX void CSYMAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZSYMAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);

_CPP_PREFIX void SSYMAMGdeleteparameters(Smat *, SILUPACKparam *, integer);
_CPP_PREFIX void DSYMAMGdeleteparameters(Dmat *, DILUPACKparam *, integer);
_CPP_PREFIX void CHERAMGdeleteparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZHERAMGdeleteparameters(Zmat *, ZILUPACKparam *, integer);
_CPP_PREFIX void CSYMAMGdeleteparameters(Cmat *, CILUPACKparam *, integer);
_CPP_PREFIX void ZSYMAMGdeleteparameters(Zmat *, ZILUPACKparam *, integer);



_CPP_PREFIX void SILUPACKparamdelete(SILUPACKparam *);
_CPP_PREFIX void DILUPACKparamdelete(DILUPACKparam *);
_CPP_PREFIX void CILUPACKparamdelete(CILUPACKparam *);
_CPP_PREFIX void ZILUPACKparamdelete(ZILUPACKparam *);


_CPP_PREFIX void Smergematrices(Smat *C, Smat *A, Smat *B, 
				integer *p, integer *invq, integer *lenB, 
				real shift, integer *buff);
_CPP_PREFIX void Dmergematrices(Dmat *C, Dmat *A, Dmat *B, 
				integer *p, integer *invq, integer *lenB, 
				doubleprecision shift, integer *buff);
_CPP_PREFIX void Cmergematrices(Cmat *C, Cmat *A, Cmat *B, 
				integer *p, integer *invq, integer *lenB, 
				complex shift, integer *buff);
_CPP_PREFIX void Zmergematrices(Zmat *C, Zmat *A, Zmat *B, 
				integer *p, integer *invq, integer *lenB, 
				doublecomplex shift, integer *buff);

_CPP_PREFIX void SSYMmergematrices(Smat *C, Smat *A, Smat *B, 
				   real *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   real shift, integer *buff);
_CPP_PREFIX void DSYMmergematrices(Dmat *C, Dmat *A, Dmat *B, 
				   doubleprecision *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   doubleprecision shift, integer *buff);
_CPP_PREFIX void CSYMmergematrices(Cmat *C, Cmat *A, Cmat *B, 
				   complex *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   complex shift, integer *buff);
_CPP_PREFIX void ZSYMmergematrices(Zmat *C, Zmat *A, Zmat *B, 
				   doublecomplex *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   doublecomplex shift, integer *buff);
_CPP_PREFIX void CHERmergematrices(Cmat *C, Cmat *A, Cmat *B, 
				   complex *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   complex shift, integer *buff);
_CPP_PREFIX void ZHERmergematrices(Zmat *C, Zmat *A, Zmat *B, 
				   doublecomplex *rowscale,
				   integer *p, integer *invq, integer *lenB, 
				   doublecomplex shift, integer *buff);

_CPP_PREFIX void SSYMAMGextractblockdiag(Smat, SAMGlevelmat *, SILUPACKparam *);
_CPP_PREFIX void DSYMAMGextractblockdiag(Dmat, DAMGlevelmat *, DILUPACKparam *);
_CPP_PREFIX void CSYMAMGextractblockdiag(Cmat, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZSYMAMGextractblockdiag(Zmat, ZAMGlevelmat *, ZILUPACKparam *);
_CPP_PREFIX void CHERAMGextractblockdiag(Cmat, CAMGlevelmat *, CILUPACKparam *);
_CPP_PREFIX void ZHERAMGextractblockdiag(Zmat, ZAMGlevelmat *, ZILUPACKparam *);


_CPP_PREFIX void DSYMdecoupleconstraints(Dmat A, integer *ia, integer **ja, 
					 doubleprecision **a, DAMGlevelmat *current, 
					 integer *indicator, doubleprecision *testvector,
					 integer *ibuff, doubleprecision *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void SSYMdecoupleconstraints(Smat A, integer *ia, integer **ja, real **a, 
					 SAMGlevelmat *current, integer *indicator, 
					 real *testvector, integer *ibuff, real *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);
_CPP_PREFIX void ZSYMdecoupleconstraints(Zmat A, integer *ia, integer **ja, 
					 doublecomplex **a, ZAMGlevelmat *current, 
					 integer *indicator, doublecomplex *testvector, 
					 integer *ibuff, doublecomplex *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void CSYMdecoupleconstraints(Cmat A, integer *ia, integer **ja, 
					 complex **a, CAMGlevelmat *current, 
					 integer *indicator, complex *testvector, 
					 integer *ibuff, complex *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);
_CPP_PREFIX void ZHERdecoupleconstraints(Zmat A, integer *ia, integer **ja, 
					 doublecomplex **a, ZAMGlevelmat *current, 
					 integer *indicator, doublecomplex *testvector, 
					 integer *ibuff, doublecomplex *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void CHERdecoupleconstraints(Cmat A, integer *ia, integer **ja, 
					 complex **a, CAMGlevelmat *current, 
					 integer *indicator, complex *testvector, 
					 integer *ibuff, complex *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);



_CPP_PREFIX void DSYMdecoupleconstraintshh(Dmat A, integer *ia, integer **ja, 
					 doubleprecision **a, DAMGlevelmat *current, 
					 integer *indicator, doubleprecision *testvector,
					 integer *ibuff, doubleprecision *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void SSYMdecoupleconstraintshh(Smat A, integer *ia, integer **ja, real **a, 
					 SAMGlevelmat *current, integer *indicator, 
					 real *testvector, integer *ibuff, real *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);
_CPP_PREFIX void ZSYMdecoupleconstraintshh(Zmat A, integer *ia, integer **ja, 
					 doublecomplex **a, ZAMGlevelmat *current, 
					 integer *indicator, doublecomplex *testvector, 
					 integer *ibuff, doublecomplex *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void CSYMdecoupleconstraintshh(Cmat A, integer *ia, integer **ja, 
					 complex **a, CAMGlevelmat *current, 
					 integer *indicator, complex *testvector, 
					 integer *ibuff, complex *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);
_CPP_PREFIX void ZHERdecoupleconstraintshh(Zmat A, integer *ia, integer **ja, 
					 doublecomplex **a, ZAMGlevelmat *current, 
					 integer *indicator, doublecomplex *testvector, 
					 integer *ibuff, doublecomplex *dbuff, 
					 doubleprecision condest, doubleprecision droptolc,
					 doubleprecision droptolS, integer PILUCparam);
_CPP_PREFIX void CHERdecoupleconstraintshh(Cmat A, integer *ia, integer **ja, 
					 complex **a, CAMGlevelmat *current, 
					 integer *indicator, complex *testvector, 
					 integer *ibuff, complex *dbuff, 
					 real condest, real droptolc, real droptolS, 
					 integer PILUCparam);


_CPP_PREFIX void silupackcast(SILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void dilupackcast(DILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void cilupackcast(CILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void zilupackcast(ZILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void silupackcastback(SILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void dilupackcastback(DILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void cilupackcastback(CILUPACKparam *options, ILUPACKparam *param);
_CPP_PREFIX void zilupackcastback(ZILUPACKparam *options, ILUPACKparam *param);



_CPP_PREFIX void write_globals(double *, size_t *);
_CPP_PREFIX void read_globals(double *, size_t *);




/* global variable to measure the timings of the components within ILUPACK */
_CPP_PREFIX void ilupack_dummy_000();
#define ILUPACK_secnds_length   12
#define ILUPACK_mem_length      12
#define ILUPACK_max_threads    256

#ifdef  _DECLARE_ILUPACK_GLOBALS_
    double ILUPACK_secnds[ILUPACK_max_threads][ILUPACK_secnds_length];
    double ILUPACK_sum_secnds[ILUPACK_secnds_length];

    size_t ILUPACK_mem[ILUPACK_max_threads][ILUPACK_mem_length];
    size_t ILUPACK_sum_mem[ILUPACK_mem_length];

/*
    ILUPACK_mem[threadnum][0], ILUPACK_sum_mem[0]:
       maximum amount of memory of auxiliary buffer "ibuff"

    ILUPACK_mem[threadnum][1], ILUPACK_sum_mem[1]:
       maximum amount of memory of auxiliary buffer "dbuff"

    ILUPACK_mem[threadnum][2], ILUPACK_sum_mem[2]:
       current amount of memory of array "jlu"
       This is used to carry all levels of the preconditioner
       inside one long integer array. Particulary, jlu consists of
       "pointers", i.e., integer references to columns of the ILU factor
       L as well as the nonzero index information

    ILUPACK_mem[threadnum][3], ILUPACK_sum_mem[3]:
       current amount of memory of array "alu"
       This is used to carry all levels of the preconditioner
       inside one long integer array. Particulary, alu consists of
       numerical values, i.e., floating point numbers of the ILU factors
       D, L (and possibly U), where the entries of D are in front while
       the entries of L refere to the nonzero index structure of jlu
       at the same positions.

    ILUPACK_mem[threadnum][4], ILUPACK_sum_mem[4]:
       peak (i.e. highest intermediate) amount of memory of array "jlu"

    ILUPACK_mem[threadnum][5], ILUPACK_sum_mem[5]:
       peak (i.e. highest intermediate) amount of memory of array "alu"

    ILUPACK_mem[threadnum][6], ILUPACK_sum_mem[6]:
       current accumulated number of pointer memory (number of rows +1)
       of the off-diagonal part F (lower triangular part)

    ILUPACK_mem[threadnum][7], ILUPACK_sum_mem[7]:
       current accumulated number nonzeros of the off-diagonal part F 
       (lower triangular part)

    ILUPACK_mem[threadnum][8], ILUPACK_sum_mem[8]:
       current accumulated memory for permutation vectors

    ILUPACK_mem[threadnum][9], ILUPACK_sum_mem[9]:
       current accumulated memory for diagonal scalings (left and right if 
       necessary)

    ILUPACK_mem[threadnum][10], ILUPACK_sum_mem[10]:
       current accumulated memory for storing integer data of coarse grid 
       systems

    ILUPACK_mem[threadnum][11], ILUPACK_sum_mem[11]:
       current accumulated memory for storing floating point data of coarse
       grid systems
 */

#else
    extern double ILUPACK_secnds[ILUPACK_max_threads][ILUPACK_secnds_length];
    extern size_t ILUPACK_mem[ILUPACK_max_threads][ILUPACK_mem_length];
    extern double ILUPACK_sum_secnds[ILUPACK_secnds_length];
    extern size_t ILUPACK_sum_mem[ILUPACK_mem_length];
#endif


void Dbisinit(integer *ipar, double *fpar, integer *wksize, integer *dsc,
	      integer *lp, integer *rp, double *wk);
void Zbisinit(integer *ipar, double *fpar, integer *wksize, integer *dsc,
	      integer *lp, integer *rp, doublecomplex *wk);

void Dmgsro(integer *full, integer *lda, integer *n, integer *m, integer *ind,
	    double *ops, double *vec, double *hh, integer *ierr);
void Zmgsro(integer *full, integer *lda, integer *n, integer *m, integer *ind,
	    double *ops, doublecomplex *vec, doublecomplex *hh, integer *ierr);

integer Dbrkdn(doubleprecision *alpha, integer *ipar);
integer Zbrkdn(doubleprecision *alpha, integer *ipar);

integer Dstopbis(integer *n,integer *ipar,integer *mvpi,
		 doubleprecision *fpar,doubleprecision *r,
		 doubleprecision *delx,doubleprecision *sx);
integer Zstopbis(integer *n,integer *ipar,integer *mvpi,
		 doubleprecision *fpar,doublecomplex *r,
		 doublecomplex *delx,doubleprecision *sx);
void Dtidycg(integer *n,integer *ipar,doubleprecision *fpar,
	     doubleprecision *sol,doubleprecision *delx);
void Ztidycg(integer *n,integer *ipar,doubleprecision *fpar,
	     doublecomplex *sol,  doublecomplex *delx);

#endif /* _ILU_PACK_H */





