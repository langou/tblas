#include "blas.h"
#include "nrm2.h"

using tblas::nrm2;

double dznrm2_(const int &n, complex<double> *x, const int &incx)
{
    if((n>0)&&(incx>0))
        return nrm2(n,x,incx);
    else
        return 0.0;
}
