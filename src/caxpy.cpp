#include "blas.h"
#include "axpy.h"

using tblas::axpy;

void caxpy_(const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    const complex<float> zero(0.0f);
    if((n>0)&&(alpha!=zero))
        axpy(n,alpha,x,incx,y,incy);
}
