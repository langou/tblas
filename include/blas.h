//
//  blas.h
//
//  Purpose
//  =======
//
//  Definition for gfortran-like legacy BLAS interface.
//

#ifndef __blas__
#define __blas__

#include <complex>

using std::complex;

extern "C"
{
    void xerbla_(const char *name, const int &info);
    
    void srotm_(const int &n, float *x, const int &incx, float *y, const int &incy, float *param);
    void srotmg_(float &d1, float &d2, float &x1, float &y1, float *param);
    void drotm_(const int &n, double *x, const int &incx, double *y, const int &incy, double *param);
    void drotmg_(double &d1, double &d2, double &x1, double &y1, double *param);
    
    void saxpy_(const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy);
    void daxpy_(const int &n, const double &alpha , double *x, const int &incx, double *y, const int &incy);
    void caxpy_(const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    void zaxpy_(const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
    
    int isamax_(const int &n, float *x, const int &incx);
    int idamax_(const int &n, double *x, const int &incx);
    int icamax_(const int &n, complex<float> *x, const int &incx);
    int izamax_(const int &n, complex<double> *x, const int &incx);
    
    float sasum_(const int &n, float *x, const int &incx);
    double dasum_(const int &n, double *x, const int &incx);
    float scasum_(const int &n, complex<float> *x, const int &incx);
    double dzasum_(const int &n, complex<double> *x, const int &incx);
    
    void scopy_(const int &n, float *x, const int &incx, float *y, const int &incy);
    void dcopy_(const int &n, double *x, const int &incx, double *y, const int &incy);
    void ccopy_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    void zcopy_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
    
    float snrm2_(const int &n, float *x, const int &incx);
    double dnrm2_(const int &n, double *x, const int &incx);
    float scnrm2_(const int &n, complex<float> *x, const int &incx);
    double dznrm2_(const int &n, complex<double> *x, const int &incx);
    
    void srot_(const int &n, float *x, const int &incx, float *y, const int &incy, const float &c, const float &s);
    void drot_(const int &n, double *x, const int &incx, double *y, const int &incy, const double &c, const double &s);
    void csrot_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy, const float &c, const float &s);
    void zdrot_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy, const double &c, const double &s);
    
    void srotg_(float &a, float &b, float &c, float &s);
    void drotg_(double &a, double &b, double &c, double &s);
    void crotg_(complex<float> &a,const complex<float> &b, float &c, complex<float> &s);
    void zrotg_(complex<double> &a,const complex<double> &b, double &c, complex<double> &s);
    
    void sscal_(const int &n, const float &a, float *x, const int &incx);
    void dscal_(const int &n, const double &a, double *x, const int &incx);
    void cscal_(const int &n, const complex<float> &a, complex<float> *x, const int &incx);
    void csscal_(const int &n, const float &a, complex<float> *x, const int &incx);
    void zscal_(const int &n, const complex<double> &a, complex<double> *x, const int &incx);
    void zdscal_(const int &n, const double &a, complex<double> *x, const int &incx);
    
    void sswap_(const int &n, float *x, const int &incx, float *y, const int &incy);
    void dswap_(const int &n, double *x, const int &incx, double *y, const int &incy);
    void cswap_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    void zswap_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
    
    float sdot_(const int &n, float *x, const int &incx, float *y, const int &incy);
    double ddot_(const int &n, double *x, const int &incx, double *y, const int &incy);
    float sdsdot_(const int &n, const float &B, float *x, const int &incx, float *y, const int &incy);
    double dsdot_(const int &n, float *x, const int &incx, float *y, const int &incy);
#ifdef __INTEL_COMPILER
    void cdotc_(complex<float> &dot, const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    void cdotu_(complex<float> &dot, const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    void zdotc_(complex<double> &dot, const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
    void zdotu_(complex<double> &dot, const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
#else
    complex<float> cdotc_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    complex<float> cdotu_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy);
    complex<double> zdotc_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
    complex<double> zdotu_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy);
#endif
    void sgemv_(const char &trans, const int &m, const int &n, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy);
    void dgemv_(const char &trans, const int &m, const int &n, const double &alpha, double *A, const int &ldA, double *x, const int &incx, const double &beta, double *y, const int &incy);
    void cgemv_(const char &trans, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy);
    void zgemv_(const char &trans, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy);
    
    void sgbmv_(const char &trans, const int &m, const int &n, const int &KL, const int &KU, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy);
    void dgbmv_(const char &trans, const int &m, const int &n, const int &KL, const int &KU, const double &alpha, double *A, const int &ldA, double *x, const int &incx, const double &beta, double *y, const int &incy);
    void cgbmv_(const char &trans, const int &m, const int &n, const int &KL, const int &KU, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy);
    void zgbmv_(const char &trans, const int &m, const int &n, const int &KL, const int &KU, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy);
    
    void ssymv_(const char &uplo, const int &n, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy);
    void dsymv_(const char &uplo, const int &n, const double &alpha, double *A, const int &ldA, double *x, const int &incx, const double &beta, double *y, const int &incy);
    void chemv_(const char &uplo, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy);
    void zhemv_(const char &uplo, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy);
    
    void ssbmv_(const char &uplo, const int &n, const int &K, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy);
    void dsbmv_(const char &uplo, const int &n, const int &K, const double &alpha, double *A, const int &ldA, double *x, const int &incx, const double &beta, double *y, const int &incy);
    void chbmv_(const char &uplo, const int &n, const int &K, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy);
    void zhbmv_(const char &uplo, const int &n, const int &K, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy);
    
    void sspmv_(const char &uplo, const int &n, const float &alpha, float *A, float *x, const int &incx, const float &beta, float *y, const int &incy);
    void dspmv_(const char &uplo, const int &n, const double &alpha, double *A, double *x, const int &incx, const double &beta, double *y, const int &incy);
    void chpmv_(const char &uplo, const int &n, const complex<float> &alpha, complex<float> *A, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy);
    void zhpmv_(const char &uplo, const int &n, const complex<double> &alpha, complex<double> *A, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy);
    
    void strmv_(const char &uplo, const char &trans, const char &diag, const int &n, float *A, const int &ldA, float *x, const int &incx);
    void dtrmv_(const char &uplo, const char &trans, const char &diag, const int &n, double *A, const int &ldA, double *x, const int &incx);
    void ctrmv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<float> *A, const int &ldA, complex<float> *x, const int &incx);
    void ztrmv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<double> *A, const int &ldA, complex<double> *x, const int &incx);
    
    void stbmv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, float *A, const int &ldA, float *x, const int &incx);
    void dtbmv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, double *A, const int &ldA, double *x, const int &incx);
    void ctbmv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, complex<float> *A, const int &ldA, complex<float> *x, const int &incx);
    void ztbmv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, complex<double> *A, const int &ldA, complex<double> *x, const int &incx);
    
    void stpmv_(const char &uplo, const char &trans, const char &diag, const int &n, float *A, float *x, const int &incx);
    void dtpmv_(const char &uplo, const char &trans, const char &diag, const int &n, double *A, double *x, const int &incx);
    void ctpmv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<float> *A, complex<float> *x, const int &incx);
    void ztpmv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<double> *A, complex<double> *x, const int &incx);
    
    void strsv_(const char &uplo, const char &trans, const char &diag, const int &n, float *A, const int &ldA, float *x, const int &incx);
    void dtrsv_(const char &uplo, const char &trans, const char &diag, const int &n, double *A, const int &ldA, double *x, const int &incx);
    void ctrsv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<float>  *A, const int &ldA, complex<float> *x, const int &incx);
    void ztrsv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<double>  *A, const int &ldA, complex<double> *x, const int &incx);
    
    void stbsv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, float *A, const int &ldA, float *x, const int &incx);
    void dtbsv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, double *A, const int &ldA, double *x, const int &incx);
    void ctbsv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, complex<float> *A, const int &ldA, complex<float> *x, const int &incx);
    void ztbsv_(const char &uplo, const char &trans, const char &diag, const int &n, const int &K, complex<double> *A, const int &ldA, complex<double> *x, const int &incx);
    
    void stpsv_(const char &uplo, const char &trans, const char &diag, const int &n, float *A, float *x, const int &incx);
    void dtpsv_(const char &uplo, const char &trans, const char &diag, const int &n, double *A, double *x, const int &incx);
    void ctpsv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<float> *A, complex<float> *x, const int &incx);
    void ztpsv_(const char &uplo, const char &trans, const char &diag, const int &n, complex<double> *A, complex<double> *x, const int &incx);
    
    void sger_(const int &m, const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy, float *A, const int &ldA);
    void dger_(const int &m, const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy, double *A, const int &ldA);
    void cgerc_(const int &m, const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy, complex<float> *A, const int &ldA);
    void zgerc_(const int &m, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A, const int &ldA);
    void cgeru_(const int &m, const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy, complex<float> *A, const int &ldA);
    void zgeru_(const int &m, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A, const int &ldA);
    
    void ssyr_(const char &uplo, const int &n, const float &alpha, float *x, const int &incx, float *A, const int &ldA);
    void dsyr_(const char &uplo, const int &n, const double &alpha, double *x, const int &incx, double *A, const int &ldA);
    void cher_(const char &uplo, const int &n, const float &alpha, complex<float> *x, const int &incx, complex<float> *A, const int &ldA);
    void zher_(const char &uplo, const int &n, const double &alpha, complex<double> *x, const int &incx, complex<double> *A, const int &ldA);
    
    void sspr_(const char &uplo, const int &n, const float &alpha, float *x, const int &incx, float *A);
    void dspr_(const char &uplo, const int &n, const double &alpha, double *x, const int &incx, double *A);
    void chpr_(const char &uplo, const int &n, const float &alpha, complex<float> *x, const int &incx, complex<float> *A);
    void zhpr_(const char &uplo, const int &n, const double &alpha, complex<double> *x, const int &incx, complex<double> *A);
    
    void ssyr2_(const char &uplo, const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy, float *A, const int &ldA);
    void dsyr2_(const char &uplo, const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy, double *A, const int &ldA);
    void cher2_(const char &uplo, const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy, complex<float> *A, const int &ldA);
    void zher2_(const char &uplo, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A, const int &ldA);
    
    void sspr2_(const char &uplo, const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy, float *A);
    void dspr2_(const char &uplo, const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy, double *A);
    void chpr2_(const char &uplo, const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy, complex<float> *A);
    void zhpr2_(const char &uplo, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A);
    
    void sgemm_(const char &transA, const char &transB, const int &m, const int &n, const int &K, const float &alpha, float *A, const int &ldA, float *B, const int &ldB, const float &beta, float *C, const int &ldC);
    void dgemm_(const char &transA, const char &transB, const int &m, const int &n, const int &K, const double &alpha, double *A, const int &ldA, double *B, const int &ldB, const double &beta, double *C, const int &ldC);
    void cgemm_(const char &transA, const char &transB, const int &m, const int &n, const int &K, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const complex<float> &beta, complex<float> *C, const int &ldC);
    void chemm_(const char &side, const char &uplo, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const complex<float> &beta, complex<float> *C, const int &ldC);
    void zgemm_(const char &transA, const char &transB, const int &m, const int &n, const int &K, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC);
    void zhemm_(const char &side, const char &uplo, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC);
    
    void ssymm_(const char &side, const char &uplo, const int &m, const int &n, const float &alpha, float *A, const int &ldA, float *B, const int &ldB, const float &beta, float *C, const int &ldC);
    void dsymm_(const char &side, const char &uplo, const int &m, const int &n, const double &alpha, double *A, const int &ldA, double *B, const int &ldB, const double &beta, double *C, const int &ldC);
    void csymm_(const char &side, const char &uplo, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const complex<float> &beta, complex<float> *C, const int &ldC);
    void zsymm_(const char &side, const char &uplo, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC);
    
    void ssyrk_(const char &uplo, const char &trans, const int &n, const int &K, const float &alpha, float *A, const int &ldA, const float &beta, float *C, const int &ldC);
    void dsyrk_(const char &uplo, const char &trans, const int &n, const int &K, const double &alpha, double *A, const int &ldA, const double &beta, double *C, const int &ldC);
    void csyrk_(const char &uplo, const char &trans, const int &n, const int &K, const complex<float> &alpha, complex<float> *A, const int &ldA, const complex<float> &beta, complex<float> *C, const int &ldC);
    void cherk_(const char &uplo, const char &trans, const int &n, const int &K, const float &alpha, complex<float> *A, const int &ldA, const float &beta, complex<float> *C, const int &ldC);
    void zsyrk_(const char &uplo, const char &trans, const int &n, const int &K, const complex<double> &alpha, complex<double> *A, const int &ldA, const complex<double> &beta, complex<double> *C, const int &ldC);
    void zherk_(const char &uplo, const char &trans, const int &n, const int &K, const double &alpha, complex<double> *A, const int &ldA, const double &beta, complex<double> *C, const int &ldC);
    
    void ssyr2k_(const char &uplo, const char &trans, const int &n, const int &K, const float &alpha, float *A, const int &ldA, float *B, const int &ldB, const float &beta, float *C, const int &ldC);
    void dsyr2k_(const char &uplo, const char &trans, const int &n, const int &K, const double &alpha, double *A, const int &ldA, double *B, const int &ldB, const double &beta, double *C, const int &ldC);
    void csyr2k_(const char &uplo, const char &trans, const int &n, const int &K, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const complex<float> &beta, complex<float> *C, const int &ldC);
    void cher2k_(const char &uplo, const char &trans, const int &n, const int &K, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const float &beta, complex<float> *C, const int &ldC);
    void zsyr2k_(const char &uplo, const char &trans, const int &n, const int &K, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC);
    void zher2k_(const char &uplo, const char &trans, const int &n, const int &K, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const double &beta, complex<double> *C, const int &ldC);
    
    void strmm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const float &alpha, float *A, const int &ldA, float *B, const int &ldB);
    void dtrmm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const double &alpha, double *A, const int &ldA, double *B, const int &ldB);
    void ctrmm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB);
    void ztrmm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB);
    
    void strsm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const float &alpha, float *A, const int &ldA, float *B, const int &ldB);
    void dtrsm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const double &alpha, double *A, const int &ldA, double *B, const int &ldB);
    void ctrsm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB);
    void ztrsm_(const char &side, const char &uplo, const char &trans, const char &diag, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB);
}
#endif
