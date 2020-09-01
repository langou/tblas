//
//  tpmv.h
//
//  Purpose
//  =======
//
//  Performs the matrix-vector product
//
//      x <- A   * x  [trans='N']
//
//      x <- A^T * x  [trans='T']
//
//      x <- A^H * x  [trans='C']
//
//  where alpha and beta are scalars, x and y are vectors and A is an
//  n-by-n triangular matrix in packaged storage format.
//
//  Arguments
//  ==========
//
//  uplo    specifies whether A is upper ('U') or lower ('L') triangular
//
//  trans   specifies the transpose operation for A as above
//
//  diag    specifies whether the matrix A is unit triangular ('U') or not ('N')
//
//  n       specifies the number of columns of the matrix A
//
//  A       triangular matrix in packaged format stored as a vector of length n*(n+1)/2
//
//  x       vector of length n
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//

#ifndef __tpmv__
#define __tpmv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void tpmv(char uplo, char trans, char diag, size_t n, T *A, T *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(n>0)
        {
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
                            A+=j+1;
                        }
                    }
                    else if(uplo=='L')
                    {
                        A+=n*(n+1)/2;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=n-j;
                            for(ptrdiff_t i=n-1;i>j;i--)
                                x[i]+=x[j]*A[i-j];
                            if(nounit)
                                x[j]*=A[0];
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[j]*=A[0];
                            for(size_t i=j+1;i<n;i++)
                                x[j]+=A[i-j]*x[i];
                            A+=n-j;
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
                            A+=j+1;
                            jx+=incx;
                        }
                    }
                    else if(uplo=='L')
                    {
                        A+=n*(n+1)/2;
                        kx+=n*incx;
                        size_t jx=kx;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=n-j;
                            jx-=incx;
                            size_t ix=kx;
                            for(ptrdiff_t i=n-1;i>j;i--)
                            {
                                ix-=incx;
                                x[ix]+=x[jx]*A[i-j];
                            }
                            if(nounit)
                                x[jx]*=A[0];
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        size_t jx=kx+n*incx;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[jx]*=A[0];
                            for(size_t i=j+1;i<n;i++)
                            {
                                ix+=incx;
                                x[jx]+=A[i-j]*x[ix];
                            }
                            A+=n-j;
                            jx+=incx;
                        }
                    }
                }
            }
        }
    }
    
    template <typename T>
    void tpmv(char uplo, char trans, char diag, size_t n, complex<T> *A, complex<T> *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(n>0)
        {
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
                            A+=j+1;
                        }
                    }
                    else if(uplo=='L')
                    {
                        A+=n*(n+1)/2;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=n-j;
                            for(ptrdiff_t i=n-1;i>j;i--)
                                x[i]+=x[j]*A[i-j];
                            if(nounit)
                                x[j]*=A[0];
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[j]*=A[0];
                            for(size_t i=j+1;i<n;i++)
                                x[j]+=A[i-j]*x[i];
                            A+=n-j;
                        }
                    }
                }
                else if(trans=='C')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[j]*=conj(A[0]);
                            for(size_t i=j+1;i<n;i++)
                                x[j]+=conj(A[i-j])*x[i];
                            A+=n-j;
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
                            A+=j+1;
                            jx+=incx;
                        }
                    }
                    else if(uplo=='L')
                    {
                        A+=n*(n+1)/2;
                        kx+=n*incx;
                        size_t jx=kx;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=n-j;
                            jx-=incx;
                            size_t ix=kx;
                            for(ptrdiff_t i=n-1;i>j;i--)
                            {
                                ix-=incx;
                                x[ix]+=x[jx]*A[i-j];
                            }
                            if(nounit)
                                x[jx]*=A[0];
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        size_t jx=kx+n*incx;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[jx]*=A[0];
                            for(size_t i=j+1;i<n;i++)
                            {
                                ix+=incx;
                                x[jx]+=A[i-j]*x[ix];
                            }
                            A+=n-j;
                            jx+=incx;
                        }
                    }
                }
                else if(trans=='C')
                {
                    if(uplo=='U')
                    {
                        A+=n*(n+1)/2;
                        size_t jx=kx+n*incx;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            A-=j+1;
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
                                x[jx]*=conj(A[0]);
                            for(size_t i=j+1;i<n;i++)
                            {
                                ix+=incx;
                                x[jx]+=conj(A[i-j])*x[ix];
                            }
                            A+=n-j;
                            jx+=incx;
                        }
                    }
                }
            }
        }
    }
}
#endif
