#include "blas.h"
#include "swap.h"

using tblas::swap;

void zswap_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy)
{
    if(n>0)
        swap(n,x,incx,y,incy);
}
