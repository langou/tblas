#include "blas.h"
#include "axpy.h"

using tblas::axpy;

void saxpy_(const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy)
{
    const float zero(0.0f);
    if((n>0)&&(alpha!=zero))
        axpy(n,alpha,x,incx,y,incy);
}
