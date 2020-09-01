#include "blas.h"
#include "axpy.h"

using tblas::axpy;

void daxpy_(const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy)
{
    const double zero(0.0);
    if((n>0)&&(alpha!=zero))
        axpy(n,alpha,x,incx,y,incy);
}
