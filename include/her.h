//
//  her.h
//
//  Purpose
//  =======
//
//  Performs a vector outer product of a vector x with its conjugate transpose
//
//      A <- alpha * x * x^H + A
//
//  where A is a Hermitian matrix, x and y are vectors, and alpha is a scalar.
//
//  Arguments
//  =========
//
//  uplo    specifies whether A is accessed as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the Hermitian matrix A
//
//  alpha   scalar multiple of the vector outer project
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  A       Hermitian matrix of order n
//
//  ldA     specifies the column length of A, must be at least n
//

#ifndef __her__
#define __her__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void her(char uplo, size_t n, T alpha, complex<T> *x, ptrdiff_t incx, complex<T> *A, size_t ldA)
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
                    A+=ldA;
                }
            }
            else if(uplo=='L')
            {
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[j]);
                    A[j]=real(A[j])+real(t*x[j]);
                    for(size_t i=j+1;i<n;i++)
                        A[i]+=x[i]*t;
                    A+=ldA;
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
                    A+=ldA;
                    jx+=incx;
                }
            }
            else if(uplo=='L')
            {
                size_t jx=kx;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*conj(x[jx]);
                    A[j]=real(A[j])+real(t*x[jx]);
                    size_t ix=jx;
                    for(size_t i=j+1;i<n;i++)
                    {
                        ix+=incx;
                        A[i]+=x[ix]*t;
                    }
                    jx+=incx;
                    A+=ldA;
                }
            }
        }
    }
}
#endif

