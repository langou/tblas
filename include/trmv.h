//
//  trmv.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-vector operations:
//
//      x <- A   * x  [trans='N']
//
//      x <- A^T * x  [trans='T']
//
//      x <- A^H * x  [trans='C']
//
//  where x is a vector and A is an n-by-n band triangular matrix.
//  If diag='U', the matrix is assumed to be unit triangular,
//  or diag='N' for nonunit triangular.
//
//  Arguments
//  ==========
//
//  uplo    specifies whether matrix is upper ('U') or lower ('L') triangular
//
//  trans   specifies the transpose operation for A as above
//
//  diag    specifies whether the matrix A is unit triangular ('U') or not ('N')
//
//  n       specifies the order of the triangular matrix A
//
//  A       triangular matrix stored as ldA-by-n array
//
//  ldA     specifies the column length of A, must be at least n
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//

#ifndef __trmv__
#define __trmv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void trmv(char uplo, char trans, char diag, size_t n, T *A, size_t ldA, T *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<j;i++)
                            x[i]+=x[j]*A[i];
                        if(nounit)
                            x[j]*=A[j];
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        for(ptrdiff_t i=n-1;i>j;i--)
                            x[i]+=x[j]*A[i];
                        if(nounit)
                            x[j]*=A[j];
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        if(nounit)
                            x[j]*=A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                            x[j]+=A[i]*x[i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]*=A[j];
                        for(size_t i=j+1;i<n;i++)
                            x[j]+=A[i]*x[i];
                        A+=ldA;
                    }
                }
            }
        }
        else
        {
            size_t kx=(incx>0)?0:(1-n)*incx;
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[ix]+=x[jx]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]*=A[j];
                        A+=ldA;
                        jx+=incx;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    kx+=n*incx;
                    size_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        size_t ix=kx;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[ix]+=x[jx]*A[i];
                        }
                        if(nounit)
                            x[jx]*=A[j];
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    size_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                        {
                            ix-=incx;
                            x[jx]+=A[i]*x[ix];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=A[j];
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[jx]+=A[i]*x[ix];
                        }
                        A+=ldA;
                        jx+=incx;
                    }
                }
            }
        }
    }

    template <typename T>
    void trmv(char uplo, char trans, char diag, size_t n, complex<T> *A, size_t ldA, complex<T> *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<j;i++)
                            x[i]+=x[j]*A[i];
                        if(nounit)
                            x[j]*=A[j];
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        for(ptrdiff_t i=n-1;i>j;i--)
                            x[i]+=x[j]*A[i];
                        if(nounit)
                            x[j]*=A[j];
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        if(nounit)
                            x[j]*=A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                            x[j]+=A[i]*x[i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]*=A[j];
                        for(size_t i=j+1;i<n;i++)
                            x[j]+=A[i]*x[i];
                        A+=ldA;
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        if(nounit)
                            x[j]*=conj(A[j]);
                        for(ptrdiff_t i=j-1;i>=0;i--)
                            x[j]+=conj(A[i])*x[i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]*=conj(A[j]);
                        for(size_t i=j+1;i<n;i++)
                            x[j]+=conj(A[i])*x[i];
                        A+=ldA;
                    }
                }
            }
        }
        else
        {
            size_t kx=(incx>0)?0:(1-n)*incx;
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[ix]+=x[jx]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]*=A[j];
                        A+=ldA;
                        jx+=incx;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    kx+=n*incx;
                    size_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        size_t ix=kx;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[ix]+=x[jx]*A[i];
                        }
                        if(nounit)
                            x[jx]*=A[j];
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    size_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                        {
                            ix-=incx;
                            x[jx]+=A[i]*x[ix];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=A[j];
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[jx]+=A[i]*x[ix];
                        }
                        A+=ldA;
                        jx+=incx;
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    size_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=conj(A[j]);
                        for(ptrdiff_t i=j-1;i>=0;i--)
                        {
                            ix-=incx;
                            x[jx]+=conj(A[i])*x[ix];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=jx;
                        if(nounit)
                            x[jx]*=conj(A[j]);
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[jx]+=conj(A[i])*x[ix];
                        }
                        A+=ldA;
                        jx+=incx;
                    }
                }
            }
        }
    }
}
#endif
