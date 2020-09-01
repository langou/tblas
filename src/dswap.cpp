#include "blas.h"
#include "swap.h"

using tblas::swap;

void dswap_(const int &n, double *x, const int &incx, double *y, const int &incy)
{
    if(n>0)
        swap(n,x,incx,y,incy);
}
