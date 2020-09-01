#include "blas.h"
#include "copy.h"

using tblas::copy;

void ccopy_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    if(n>0)
        copy(n,x,incx,y,incy);
}
