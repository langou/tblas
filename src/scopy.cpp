#include "blas.h"
#include "copy.h"

using tblas::copy;

void scopy_(const int &n, float *x, const int &incx, float *y, const int &incy)
{
    if(n>0)
        copy(n,x,incx,y,incy);
}
