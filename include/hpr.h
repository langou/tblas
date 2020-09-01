//
//  hpr.h
//
//  Purpose
//  =======
//
//  Performs the outer product of a vector x 
//
//      A <- alpha * x * x^H + A
//
//  where A is a Hermitian matrix in packed storage format,
//  x is a vector, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  uplo    specifies whether A is stored as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the Hermitian matrix A
//
//  alpha   scalar multiple of the vector outer project
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  A       Hermitian matrix of order n stored in packed format in a
//          vector of length n(n+1)/2
//

#ifndef __hpr__
#define __hpr__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void hpr(char uplo, size_t n, T alpha, complex<T> *x, ptrdiff_t incx, complex<T> *A)
    {
        const T zero(0.0);
        
        if((n==0)||(alpha==zero))
            return;
        
        if(incx==1)
        {
            if(uplo=='U')
            {
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[j]);
                    for(size_t i=0;i<j;i++)
                        A[i]+=x[i]*t;
                    A[j]=real(A[j])+real(x[j]*t);
                    A+=j+1;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[j]);
                    A[0]=real(A[0])+real(t*x[j]);
                    for(size_t i=j+1;i<n;i++)
                        A[i-j]+=x[i]*t;
                    A+=n-j;
                }
            }
        }
        else
        {
            size_t kx=(incx>0)?0:(1-n)*incx;
            if(uplo=='U')
            {
                size_t jx=kx;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[jx]);
                    size_t ix=kx;
                    for(size_t i=0;i<j;i++)
                    {
                        A[i]+=x[ix]*t;
                        ix+=incx;
                    }
                    A[j]=real(A[j])+real(x[jx]*t);
                    A+=j+1;
                    jx+=incx;
                }
            }
            else if(uplo=='L')
            {
                size_t jx=kx;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[jx]);
                    A[0]=real(A[0])+real(t*x[jx]);
                    size_t ix=jx;
                    for(size_t i=j+1;i<n;i++)
                    {
                        ix+=incx;
                        A[i-j]+=x[ix]*t;
                    }
                    jx+=incx;
                    A+=n-j;
                }
            }
        }
    }
}
#endif
