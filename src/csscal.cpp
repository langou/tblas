#include "blas.h"
#include "scal.h"

using tblas::scal;

void csscal_(const int &n, const float &alpha, complex<float> *x, const int &incx)
{
    if((n>0)&&(incx>0))
        scal(n,alpha,x,incx);
}
