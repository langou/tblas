#include "blas.h"
#include "asum.h"

using tblas::asum;

float sasum_(const int &n, float *x, const int &incx)
{
    if((n>0)&&(incx>0))
        return asum(n,x,incx);
    else
        return 0.0f;
}
