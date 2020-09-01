#include "blas.h"
#include "scal.h"

using tblas::scal;

void cscal_(const int &n, const complex<float> &alpha, complex<float> *x, const int &incx)
{
    if((n>0)&&(incx>0))
        scal(n,alpha,x,incx);
}
