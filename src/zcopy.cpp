#include "blas.h"
#include "copy.h"

using tblas::copy;

void zcopy_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy)
{
    if(n>0)
        copy(n,x,incx,y,incy);
}
