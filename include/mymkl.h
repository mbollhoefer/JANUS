/* useful MKL extensions to sparse matrices
   $Id: mymkl.h 3203 2017-04-16 22:07:48Z bolle $ */
#ifndef _MYMKL_H_
#define _MYMKL_H
  
#define MKL_DOMAIN_ALL      0
#define MKL_DOMAIN_BLAS     1
#define MKL_DOMAIN_FFT      2
#define MKL_DOMAIN_VML      3
#define MKL_DOMAIN_PARDISO  4


int mkl_get_max_threads();
int mkl_get_dynamic();
int mkl_domain_get_max_threads(int domain);
int mkl_domain_set_num_threads(int nt, int domain);


void cblas_ssctr(integer N, real *X, integer *indx, real *Y);
void cblas_dsctr(integer N, doubleprecision *X, integer *indx, 
		 doubleprecision *Y);
void cblas_csctr(integer N, complex *X, integer *indx, complex *Y);
void cblas_zsctr(integer N, doublecomplex *X, integer *indx, doublecomplex *Y);

void cblas_sgthr(integer N, real *Y, real *X, integer *indx);
void cblas_dgthr(integer N, doubleprecision *Y, doubleprecision *X, 
		 integer *indx);
void cblas_cgthr(integer N, complex *Y, complex *X, integer *indx);
void cblas_zgthr(integer N, doublecomplex *Y, doublecomplex *X, integer *indx);

void cblas_saxpyi(integer N, real alpha, real *X, integer *indx, 
		  real *Y);
void cblas_daxpyi(integer N, doubleprecision alpha, doubleprecision *X, 
		  integer *indx, doubleprecision *Y);
void cblas_caxpyi(integer N, complex *alpha, complex *X, integer *indx, 
		  complex *Y);
void cblas_zaxpyi(integer N, doublecomplex *alpha, doublecomplex *X, 
		  integer *indx, doublecomplex *Y);

doubleprecision cblas_sdoti(integer N, real *X, integer *indx,
			    real *Y);
doubleprecision cblas_ddoti(integer N, doubleprecision *X, integer *indx,
			    doubleprecision *Y);
void   cblas_cdotui_sub(integer N, complex *X, integer *indx,
                        complex *Y, complex *dotui);
void   cblas_cdotci_sub(integer N, complex *X, integer *indx,
                        complex *Y, complex *dotci);
void   cblas_zdotui_sub(integer N, doublecomplex *X, integer *indx,
                        doublecomplex *Y, doublecomplex *dotui);
void   cblas_zdotci_sub(integer N, doublecomplex *X, integer *indx,
                        doublecomplex *Y, doublecomplex *dotci);



#endif

