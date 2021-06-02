#ifndef _ILU_PACK_HH
#define _ILU_PACK_HH

#define _CPP_PREFIX   extern "C"

extern "C" {
#include <ilupack1.h>
}

namespace JANUS {
  using ::SparseMatrix;
  using ::DSparseMatrix;
  using ::ZSparseMatrix;

  using ::SparseBlockMatrix;
  using ::DSparseBlockMatrix;
  using ::ZSparseBlockMatrix;

  using ::JanusOptions;
  using ::JanusPrec;
  using ::JanusDefaultOptions;
  using ::JanusDelete;
  using ::JanusInit;
  using ::JanusFactor;
  using ::JanusNnz;
  using ::JanusSolver;

  using ::JanusSol;
  using ::JanusSolH;
  using ::JanusSolT;

  using ::PrintMatrix;
  using ::DPrintMatrix;
  using ::ZPrintMatrix;

  using ::DSYMbfspai;
  using ::ZSYMbfspai;
  using ::ZHERbfspai;

  using ::DSYMselbinv;
  using ::ZSYMselbinv;
  using ::ZHERselbinv;

  using ::MatVec;
  using ::DMatVec;
  using ::ZMatVec;
  using ::DMatTVec;
  using ::ZMatTVec;
  using ::DMatHVec;
  using ::ZMatHVec;
  using ::DSYMMatVec;
  using ::ZSYMMatVec;
  using ::ZHERMatVec;

  using ::SparseBlockDelete;
  using ::DSparseBlockDelete;
  using ::ZSparseBlockDelete;
  
  using ::SparseBlockInit;
  using ::DSparseBlockInit;
  using ::ZSparseBlockInit;
  
  using ::SparseDelete;
  using ::DSparseDelete;
  using ::ZSparseDelete;
  
  using ::SparseInit;
  using ::DSparseInit;
  using ::ZSparseInit;
}

#endif
