//
//  trsv.h
//
//  Purpose
//  =======
//
//  Solves one of the following systems:
//
//      A   * x = b  [trans='N']
//
//      A^T * x = b  [trans='T']
//
//      A^H * x = b  [trans='C']
//
//  where x and b are vectors and A is an n-by-n triangular matrix.
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
//  A       triangular matrix stored in ldA-by-n array
//
//  ldA     column length of the matrix A, must be at least n
//
//  x       vector of length n containing the vector b on entry, and the solution x on exit
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//

#ifndef __trsv__
#define __trsv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void trsv(char uplo, char trans, char diag, size_t n, T *A, size_t ldA, T *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        if(nounit)
                            x[j]=x[j]/A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                            x[i]-=x[j]*A[i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]=x[j]/A[j];
                        for(size_t i=j+1;i<n;i++)
                            x[i]-=x[j]*A[i];
                        A+=ldA;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<j;i++)
                            x[j]-=x[i]*A[i];
                        if(nounit)
                            x[j]=x[j]/A[j];
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
                            x[j]-=x[i]*A[i];
                        if(nounit)
                            x[j]=x[j]/A[j];
                    }
                }
            }
        }
        else
        {
            ptrdiff_t kx=(incx>0)?0:(1-n)*incx;
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    ptrdiff_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        ptrdiff_t ix=jx;
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                        {
                            ix-=incx;
                            x[ix]-=x[jx]*A[i];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        ptrdiff_t ix=jx;
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i];
                        }
                        jx+=incx;
                        A+=ldA;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        ptrdiff_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        jx+=incx;
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    kx+=n*incx;
                    ptrdiff_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        ptrdiff_t ix=kx;
                        A-=ldA;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*A[i];
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                    }
                }
            }
        }
    }

    template <typename T>
    void trsv(char uplo, char trans, char diag, size_t n, complex<T> *A, size_t ldA, complex<T> *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        if(nounit)
                            x[j]=x[j]/A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                            x[i]-=x[j]*A[i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]=x[j]/A[j];
                        for(size_t i=j+1;i<n;i++)
                            x[i]-=x[j]*A[i];
                        A+=ldA;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<j;i++)
                            x[j]-=x[i]*A[i];
                        if(nounit)
                            x[j]=x[j]/A[j];
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
                            x[j]-=x[i]*A[i];
                        if(nounit)
                            x[j]=x[j]/A[j];
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<j;i++)
                            x[j]-=x[i]*conj(A[i]);
                        if(nounit)
                            x[j]=x[j]/conj(A[j]);
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
                            x[j]-=x[i]*conj(A[i]);
                        if(nounit)
                            x[j]=x[j]/conj(A[j]);
                    }
                }
            }
        }
        else
        {
            ptrdiff_t kx=(incx>0)?0:(1-n)*incx;
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*ldA;
                    ptrdiff_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        jx-=incx;
                        ptrdiff_t ix=jx;
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        for(ptrdiff_t i=j-1;i>=0;i--)
                        {
                            ix-=incx;
                            x[ix]-=x[jx]*A[i];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        ptrdiff_t ix=jx;
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i];
                        }
                        jx+=incx;
                        A+=ldA;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        ptrdiff_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        jx+=incx;
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    kx+=n*incx;
                    ptrdiff_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        ptrdiff_t ix=kx;
                        A-=ldA;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*A[i];
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        ptrdiff_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*conj(A[i]);
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/conj(A[j]);
                        jx+=incx;
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    kx+=n*incx;
                    ptrdiff_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        ptrdiff_t ix=kx;
                        A-=ldA;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*conj(A[i]);
                        }
                        if(nounit)
                            x[jx]=x[jx]/conj(A[j]);
                    }
                }
            }
        }
    }
}
#endif
