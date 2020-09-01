#include "blas.h"
#include "dot.h"

using tblas::dot;

double ddot_(const int &n, double *x, const int &incx, double *y, const int &incy)
{
    double sum(0.0);
    if(n>=0)
        sum=dot(n,sum,x,incx,y,incy);
    return sum;
}
