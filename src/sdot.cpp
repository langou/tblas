#include "blas.h"
#include "dot.h"

using tblas::dot;

float sdot_(const int &n, float *x, const int &incx, float *y, const int &incy)
{
    float sum(0.0f);
    if(n>=0)
        sum=dot(n,sum,x,incx,y,incy);
    return sum;
}
