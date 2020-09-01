#include "blas.h"
#include "dot.h"

using tblas::dot;

#ifdef __INTEL_COMPILER
void zdotu_(complex<double> &sum, const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy)
{
    const complex<double> zero(0.0,0.0);
    sum=zero;
#else
complex<double> zdotu_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy)
{
    const complex<double> zero(0.0,0.0);
    complex<double> sum(zero);
#endif
    if(n>=0)
        sum=dot(n,sum,x,incx,y,incy);
#ifndef __INTEL_COMPILER
    return sum;
#endif
}
