//
//  tpsv.h
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
//  where x and b are vectors and A is an n-by-n triangular matrix in packed storage
//  format.  If diag='U', the matrix is assumed to be unit triangular,
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
//  A       triangular matrix in packaged format stored as a vector of length n*(n+1)/2
//
//  x       vector of length n containing the vector b on entry, and the solution x on exit
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//

#ifndef __tpsv__
#define __tpsv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void tpsv(char uplo, char trans, char diag, size_t n, T *A, T *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*(n+1)/2;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=j+1;
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
                            x[j]=x[j]/A[0];
                        for(size_t i=j+1;i<n;i++)
                            x[i]-=x[j]*A[i-j];
                        A+=n-j;
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
                            x[j]-=x[i]*A[i-j];
                        if(nounit)
                            x[j]=x[j]/A[0];
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
                    A+=n*(n+1)/2;
                    size_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=j+1;
                        jx-=incx;
                        size_t ix=jx;
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
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                        size_t ix=jx;
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i-j];
                        }
                        jx+=incx;
                        A+=n-j;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        jx+=incx;
                        A+=j+1;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*(n+1)/2;
                    kx+=n*incx;
                    size_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        size_t ix=kx;
                        A-=n-j;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*A[i-j];
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                    }
                }
            }
        }
    }
    
    template <typename T>
    void tpsv(char uplo, char trans, char diag, size_t n, complex<T> *A, complex<T> *x, ptrdiff_t incx)
    {
        const bool nounit(diag=='N');
        
        if(incx==1)
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    A+=n*(n+1)/2;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=j+1;
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
                            x[j]=x[j]/A[0];
                        for(size_t i=j+1;i<n;i++)
                            x[i]-=x[j]*A[i-j];
                        A+=n-j;
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
                            x[j]-=x[i]*A[i-j];
                        if(nounit)
                            x[j]=x[j]/A[0];
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
                            x[j]-=x[i]*conj(A[i-j]);
                        if(nounit)
                            x[j]=x[j]/conj(A[0]);
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
                    A+=n*(n+1)/2;
                    size_t jx=kx+n*incx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=j+1;
                        jx-=incx;
                        size_t ix=jx;
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
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                        size_t ix=jx;
                        for(size_t i=j+1;i<n;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i-j];
                        }
                        jx+=incx;
                        A+=n-j;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[j];
                        jx+=incx;
                        A+=j+1;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*(n+1)/2;
                    kx+=n*incx;
                    size_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        size_t ix=kx;
                        A-=n-j;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*A[i-j];
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t ix=kx;
                        for(size_t i=0;i<j;i++)
                        {
                            x[jx]-=x[ix]*conj(A[i]);
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/conj(A[j]);
                        jx+=incx;
                        A+=j+1;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*(n+1)/2;
                    kx+=n*incx;
                    size_t jx=kx;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        jx-=incx;
                        size_t ix=kx;
                        A-=n-j;
                        for(ptrdiff_t i=n-1;i>j;i--)
                        {
                            ix-=incx;
                            x[jx]-=x[ix]*conj(A[i-j]);
                        }
                        if(nounit)
                            x[jx]=x[jx]/conj(A[0]);
                    }
                }
            }
        }
    }
}
#endif
