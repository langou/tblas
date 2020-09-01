#include "blas.h"
#include "dot.h"

using tblas::dot;

#ifdef __INTEL_COMPILER
void cdotu_(complex<float> &sum, const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    const complex<float> zero(0.0f,0.0f);
    sum=zero;
#else
complex<float> cdotu_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    const complex<float> zero(0.0f,0.0f);
    complex<float> sum(zero);
#endif
    if(n>=0)
        sum=dot(n,sum,x,incx,y,incy);
#ifndef __INTEL_COMPILER
    return sum;
#endif
}
