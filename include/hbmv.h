//
//  hbmv.h
//
//  Purpose
//  =======
//
//  Performs the matrix-vector operation:
//
//      y <- alpha * A * x + beta * y
//
//  where alpha and beta are scalars, x and y are vectors and A is an
//  n-by-n Hermitian band matrix, with k sub-diagonals and k super-diagonals.
//
//  Arguments
//  ==========
//
//  uplo    specifies whether A is stored as upper ('U') or lower ('L') triangular
//
//  n       specifies the order of the Hermitian matrix A
//
//  k       specifies the number of sub or super diagonals of the matrix A
//
//  alpha   scalar multiple of the matrix A
//
//  A       Hermitian band matrix stored as (k+1)-by-n array
//
//  ldA     specifies the column length of A, must be at least k+1
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

#ifndef __hbmv__
#define __hbmv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void hbmv(char uplo, size_t n, size_t k, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *x, ptrdiff_t incx, complex<T> beta, complex<T> *y, ptrdiff_t incy)
    {
        const complex<T> one(1.0f);
        const complex<T> zero(0.0f);
        
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
                        size_t i0=(j>k)?j-k:0;
                        for(size_t i=i0;i<j;i++)
                        {
                            y[i]+=temp*A[k-j+i];
                            sum+=conj(A[k-j+i])*x[i];
                        }
                        y[j]+=temp*real(A[k])+alpha*sum;
                        A+=ldA;
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
                        size_t i0=(j>k)?j-k:0;
                        for(size_t i=i0;i<j;i++)
                        {
                            y[iy]+=temp*A[k-j+i];
                            sum+=conj(A[k-j+i])*x[ix];
                            ix+=incx;
                            iy+=incy;
                        }
                        y[jy]+=temp*real(A[k])+alpha*sum;
                        A+=ldA;
                        jx+=incx;
                        jy+=incy;
                        if(j>=k)
                        {
                            kx+=incx;
                            ky+=incy;
                        }
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
                        size_t N=(j+k<n)?j+k+1:n;
                        for(size_t i=j+1;i<N;i++)
                        {
                            y[i]+=temp*A[i-j];
                            sum+=conj(A[i-j])*x[i];
                        }
                        y[j]+=alpha*sum;
                        A+=ldA;
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
                        size_t N=(j+k<n)?j+k+1:n;
                        y[jy]+=temp*real(A[0]);
                        for(size_t i=j+1;i<N;i++)
                        {
                            ix+=incx;
                            iy+=incy;
                            y[iy]+=temp*A[i-j];
                            sum+=conj(A[i-j])*x[ix];
                        }
                        y[jy]+=alpha*sum;
                        jx+=incx;
                        jy+=incy;
                        A+=ldA;
                    }
                }
            }
        }
    }
}
#endif
