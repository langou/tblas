#include "blas.h"
#include "scal.h"

using tblas::scal;

void sscal_(const int &n, const float &alpha, float *x, const int &incx)
{
    if((n>0)&&(incx>0))
        scal(n,alpha,x,incx);
}
