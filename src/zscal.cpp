#include "blas.h"
#include "scal.h"

using tblas::scal;

void zscal_(const int &n, const complex<double> &alpha, complex<double> *x, const int &incx)
{
    if((n>0)&&(incx>0))
        scal(n,alpha,x,incx);
}
