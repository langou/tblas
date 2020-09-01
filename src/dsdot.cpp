#include "blas.h"
#include "dot.h"

using tblas::dot;

double dsdot_(const int &n, float *x, const int &incx, float *y, const int &incy)
{
    double sum(0.0);
    if(n>=0)
        sum=dot(n,sum,x,incx,y,incy);
    return sum;
}
