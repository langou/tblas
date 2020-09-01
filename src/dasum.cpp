#include "blas.h"
#include "asum.h"

using tblas::asum;

double dasum_(const int &n, double *x, const int &incx)
{
    if((n>0)&&(incx>0))
        return asum(n,x,incx);
    else
        return 0.0;
}
