//
//  gbmv.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-vector operations:
//
//      y <- alpha * A   * x + beta * y  [trans='N']
//
//      y <- alpha * A^T * x + beta * y  [trans='T']
//
//      y <- alpha * A^H * x + beta * y  [trans='C']
//
//  where alpha and beta are scalars, x and y are vectors and A is an
//  m-by-n band matrix, with kl sub-diagonals and ku super-diagonals.
//
//  Arguments
//  ==========
//
//  trans   specifies the transpose operation for A as above
//
//  m       specifies the number of rows of the matrix A
//
//  n       specifies the number of columns of the matrix A
//
//  kl      specifies the number of sub-diagonals of the matrix A
//
//  ku      specifies the number of super-diagonals of the matrix A
//
//  alpha   scalar multiple of the matrix A
//
//  A       band matrix stored as (kl+ku+1)-by-n array
//
//  ldA     specifies the column length of A, must be at least kl+ku+1
//
//  x       vector of length n if trans='N', or length m otherwise
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//
//  beta    scalar multiple of the vector y
//
//  y       vector of length m if trans='N', or length n otherwise
//
//  incy    stride of vector y; if negative, y is stored in reverse order
//

#ifndef __gbmv__
#define __gbmv__

#include <cstddef>
#include <complex>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void gbmv(char trans, size_t m, size_t n, size_t kl, size_t ku, T alpha, T *A, size_t ldA, T *x, ptrdiff_t incx, T beta, T *y, ptrdiff_t incy)
    {
        const T one(1.0);
        const T zero(0.0);
        
        if((m==0)||(n==0)) return;

        size_t lenx=(trans=='N')?n:m;
        size_t leny=(trans=='N')?m:n;
        size_t kx=(incx>0)?0:(1-lenx)*incx;
        size_t ky=(incy>0)?0:(1-leny)*incy;
        
        if(beta==zero)
        {
            if(incy==1)
            {
                for(size_t i=0;i<leny;i++)
                    y[i]=zero;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<leny;i++)
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
                for(size_t i=0;i<leny;i++)
                    y[i]*=beta;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<leny;i++)
                {
                    y[iy]*=beta;
                    iy+=incy;
                }
            }
        }
        
        if(alpha!=zero)
        {
            if(trans=='N')
            {
                size_t jx=kx;
                if(incy==1)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        T temp=alpha*x[jx];
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                            y[i]+=temp*A[k+i];
                        A+=ldA;
                        jx+=incx;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        T temp=alpha*x[jx];
                        size_t iy=ky;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                        {
                            y[iy]+=temp*A[k+i];
                            iy+=incy;
                        }
                        A+=ldA;
                        jx+=incx;
                        if(j>=ku)
                            ky+=incy;
                    }
                }
            }
            else
            {
                size_t jy=ky;
                if(incx==1)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        T temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                            temp+=A[k+i]*x[i];
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        T temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        size_t ix=kx;
                        for(size_t i=i0;i<im;i++)
                        {
                            temp+=A[k+i]*x[ix];
                            ix+=incx;
                        }
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                        if(j>=ku)
                            kx+=incx;
                    }
                }
            }
        }
    }
    
    template <typename T>
    void gbmv(char trans, size_t m, size_t n, size_t kl, size_t ku, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *x, ptrdiff_t incx, complex<T> beta, complex<T> *y, ptrdiff_t incy)
    {
        const complex<T> one(1.0);
        const complex<T> zero(0.0);
        
        if((m==0)||(n==0)) return;

        size_t lenx=(trans=='N')?n:m;
        size_t leny=(trans=='N')?m:n;
        size_t kx=(incx>0)?0:(1-lenx)*incx;
        size_t ky=(incy>0)?0:(1-leny)*incy;
        
        if(beta==zero)
        {
            if(incy==1)
            {
                for(size_t i=0;i<leny;i++)
                    y[i]=zero;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<leny;i++)
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
                for(size_t i=0;i<leny;i++)
                    y[i]*=beta;
            }
            else
            {
                size_t iy=ky;
                for(size_t i=0;i<leny;i++)
                {
                    y[iy]*=beta;
                    iy+=incy;
                }
            }
        }
        
        if(alpha!=zero)
        {
            if(trans=='N')
            {
                size_t jx=kx;
                if(incy==1)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[jx];
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                            y[i]+=temp*A[k+i];
                        A+=ldA;
                        jx+=incx;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=alpha*x[jx];
                        size_t iy=ky;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                        {
                            y[iy]+=temp*A[k+i];
                            iy+=incy;
                        }
                        A+=ldA;
                        jx+=incx;
                        if(j>=ku)
                            ky+=incy;
                    }
                }
            }
            else if(trans=='T')
            {
                size_t jy=ky;
                if(incx==1)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                            temp+=A[k+i]*x[i];
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        size_t ix=kx;
                        for(size_t i=i0;i<im;i++)
                        {
                            temp+=A[k+i]*x[ix];
                            ix+=incx;
                        }
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                        if(j>=ku)
                            kx+=incx;
                    }
                }
            }
            else if(trans=='C')
            {
                size_t jy=ky;
                if(incx==1)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        for(size_t i=i0;i<im;i++)
                            temp+=conj(A[k+i])*x[i];
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> temp=zero;
                        size_t k=ku-j;
                        size_t i0=(j>ku)?j-ku:0;
                        size_t im=(j+kl<m)?j+kl+1:m;
                        size_t ix=kx;
                        for(size_t i=i0;i<im;i++)
                        {
                            temp+=conj(A[k+i])*x[ix];
                            ix+=incx;
                        }
                        y[jy]+=alpha*temp;
                        jy+=incy;
                        A+=ldA;
                        if(j>=ku)
                            kx+=incx;
                    }
                }
            }
        }
    }

}
#endif
