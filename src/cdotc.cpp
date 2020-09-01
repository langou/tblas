#include "blas.h"
#include "dotc.h"

using tblas::dotc;

#ifdef __INTEL_COMPILER
void cdotc_(complex<float> &sum, const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    const complex<float> zero(0.0f,0.0f);
    sum=zero;
#else
complex<float> cdotc_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy)
{
    const complex<float> zero(0.0f,0.0f);
    complex<float> sum(zero);
#endif
    if(n>=0)
        sum=dotc(n,sum,x,incx,y,incy);
#ifndef __INTEL_COMPILER
    return sum;
#endif
}
