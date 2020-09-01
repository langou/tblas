//
//  hpmv.h
//
//  Purpose
//  =======
//
//  Performs the matrix-vector operation:
//
//      y <- alpha * A * x + beta * y
//
//  where alpha and beta are scalars, x and y are vectors and A is an
//  n-by-n Hermitian matrix is packed storage format.
//
//  Arguments
//  ==========
//
//  uplo    specifies whether A is stored as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the Hermitian matrix A
//
//  alpha   scalar multiple of the matrix A
//
//  A       Hermitian matrix of order n stored in packed format in a
//          vector of length n(n+1)/2
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  beta    scalar multiple of the vector y
//
//  y       vector of length n
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//

#ifndef __hpmv__
#define __hpmv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void hpmv(char uplo, size_t n, complex<T> alpha, complex<T> *A, complex<T> *x, ptrdiff_t incx, complex<T> beta, complex<T> *y, ptrdiff_t incy)
    {
        const complex<T> one(1.0);
        const complex<T> zero(0.0);
        
        if(n==0)
            return;
        
        size_t kx=(incx>0)?0:(1-n)*incx;
        size_t ky=(incy>0)?0:(1-n)*incy;
        
        if(beta==zero)
        {
            if(incy==1)
            {
                for(size_t i=0;i<n;i++)
                    y[i]=zero;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<n;i++)
                {
                    y[iy]=zero;
                    iy+=incy;
                }
                
            }
        }
        else if(beta!=one)
        {
            if(incy==1)
            {
                for(size_t i=0;i<n;i++)
                    y[i]*=beta;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<n;i++)
                {
                    y[iy]*=beta;
                    iy+=incy;
                }
            }
        }
        
        if(alpha!=zero)
        {
            if(uplo=='U')
            {
                if((incx==1)&&(incy==1))
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[j];
                        complex<T> sum=zero;
                        for(size_t i=0;i<j;i++)
                        {
                            y[i]+=temp*A[i];
                            sum+=conj(A[i])*x[i];
                        }
                        y[j]+=temp*real(A[j])+alpha*sum;
                        A+=j+1;
                    }
                }
                else
                {
                    size_t jx=kx;
                    size_t jy=ky;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[jx];
                        complex<T> sum=zero;
                        size_t ix=kx;
                        size_t iy=ky;
                        for(size_t i=0;i<j;i++)
                        {
                            y[iy]+=temp*A[i];
                            sum+=conj(A[i])*x[ix];
                            ix+=incx;
                            iy+=incy;
                        }
                        y[jy]+=temp*real(A[j])+alpha*sum;
                        A+=j+1;
                        jx+=incx;
                        jy+=incy;
                    }
                }
            }
            else if(uplo=='L')
            {
                if((incx==1)&&(incy==1))
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[j];
                        complex<T> sum=zero;
                        y[j]+=temp*real(A[0]);
                        for(size_t i=j+1;i<n;i++)
                        {
                            y[i]+=temp*A[i-j];
                            sum+=conj(A[i-j])*x[i];
                        }
                        y[j]+=alpha*sum;
                        A+=n-j;
                    }
                }
                else
                {
                    size_t jx=kx;
                    size_t jy=ky;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[jx];
                        complex<T> sum=zero;
                        size_t ix=jx;
                        size_t iy=jy;
                        y[jy]+=temp*real(A[0]);
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            iy+=incy;
                            y[iy]+=temp*A[i-j];
                            sum+=conj(A[i-j])*x[ix];
                        }
                        y[jy]+=alpha*sum;
                        jx+=incx;
                        jy+=incy;
                        A+=n-j;
                    }
                }
            }
        }
    }
}
#endif
