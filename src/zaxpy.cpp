#include "blas.h"
#include "axpy.h"

using tblas::axpy;

void zaxpy_(const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy)
{
    const complex<double> zero(0.0);
    if((n>0)&&(alpha!=zero))
        axpy(n,alpha,x,incx,y,incy);
}
