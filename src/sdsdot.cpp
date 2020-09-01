#include "blas.h"
#include "dot.h"

using tblas::dot;

float sdsdot_(const int &n, const float &b,float *x, const int &incx, float *y, const int &incy)
{
    double sum(b);
    float fsum(0.0f);
    if(n>=0)
        fsum=static_cast<float>(dot(n,sum,x,incx,y,incy));
    return fsum;
}
