#include "blas.h"
#include "imax.h"

using tblas::imax;

int idamax_(const int &n, double *x, const int &incx)
{
    if((n<1)||(incx<1))
        return 0;
    else
        return imax(n,x,incx)+1;
}
