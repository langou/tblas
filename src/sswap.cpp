#include "blas.h"
#include "swap.h"

using tblas::swap;

void sswap_(const int &n, float *x, const int &incx, float *y, const int &incy)
{
    if(n>0)
        swap(n,x,incx,y,incy);
}
