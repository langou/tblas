#include "blas.h"
#include "scal.h"

using tblas::scal;

void dscal_(const int &n, const double &alpha, double *x, const int &incx)
{
    if((n>0)&&(incx>0))
        scal(n,alpha,x,incx);
}
