#include "blas.h"
#include "swap.h"

using tblas::swap;

void cswap_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    if(n>0)
        swap(n,x,incx,y,incy);
}
