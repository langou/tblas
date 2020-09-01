#include "blas.h"
#include "copy.h"

using tblas::copy;

void dcopy_(const int &n, double *x, const int &incx, double *y, const int &incy)
{
    if(n>0)
        copy(n,x,incx,y,incy);
}
