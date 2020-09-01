//
//  tbsv.h
//
//  Purpose
//  =======
//
//  Solves on the following systems:
//
//      A   * x = b  [trans='N']
//
//      A^T * x = b  [trans='T']
//
//      A^H * x = b  [trans='C']
//
//  where x and b are vectors and A is an n-by-n band triangular band matrix
//  with k+1 diagonals.  If diag='U', the matrix is assumed to be unit triangular,
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
//  k       specifies the number of sub/super-diagonals of the matrix A
//
//  A       band matrix stored as (k+1)-by-n array
//
//  ldA     specifies the column length of A, must be at least k+1
//
//  x       vector of length n containing the vector b on entry, and the solution x on exit
//
//  incx    stride of vector x; if negative, x is stored in reverse order
//

#ifndef __tbsv__
#define __tbsv__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    inline size_t maxsub(size_t j, size_t k)
    {
        return (j>k)?j-k:0;
    }
    
    inline size_t minadd(size_t j, size_t k, size_t n)
    {
        return (j+k<n)?j+k:n;
    }

    template <typename T>
    void tbsv(char uplo, char trans, char diag, size_t n, size_t k, T *A, size_t ldA, T *x, ptrdiff_t incx)
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
                            x[j]=x[j]/A[k];
                        ptrdiff_t i0=maxsub(j,k);
                        for(ptrdiff_t i=j-1;i>=i0;i--)
                            x[i]-=x[j]*A[k-j+i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]=x[j]/A[0];
                        size_t im=minadd(j,k,n-1);
                        for(size_t i=j+1;i<=im;i++)
                            x[i]-=x[j]*A[i-j];
                        A+=ldA;
                    }
                }
            }
            else if((trans=='T')||(trans=='C'))
            {
                if(uplo=='U')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        size_t i0=maxsub(j,k);
                        for(size_t i=i0;i<j;i++)
                            x[j]-=x[i]*A[k-j+i];
                        if(nounit)
                            x[j]=x[j]/A[k];
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        for(ptrdiff_t i=im;i>j;i--)
                            x[j]-=x[i]*A[i-j];
                        if(nounit)
                            x[j]=x[j]/A[0];
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
                            x[jx]=x[jx]/A[k];
                        ptrdiff_t i0=maxsub(j,k);
                        for(ptrdiff_t i=j-1;i>=i0;i--)
                        {
                            ix-=incx;
                            x[ix]-=x[jx]*A[k-j+i];
                        }
                    }
                }
                else if(uplo=='L')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                        size_t ix=jx;
                        size_t im=minadd(j,k,n-1);
                        for(size_t i=j+1;i<=im;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i-j];
                        }
                        jx+=incx;
                        A+=ldA;
                    }
                }
            }
            else if((trans=='T')||(trans=='C'))
            {
                if(uplo=='U')
                {
                    ptrdiff_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t i0=maxsub(j,k);
                        size_t ix=kx+i0*incx;
                        for(size_t i=i0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[k-j+i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[k];
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
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        ptrdiff_t ix=kx-(n-im)*incx;
                        for(ptrdiff_t i=im;i>j;i--)
                        {
                            x[jx]-=x[ix]*A[i-j];
                            ix-=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[0];
                    }
                }
            }
        }
    }
    
    template <typename T>
    void tbsv(char uplo, char trans, char diag, size_t n, size_t k, complex<T> *A, size_t ldA, complex<T> *x, ptrdiff_t incx)
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
                            x[j]=x[j]/A[k];
                        ptrdiff_t i0=maxsub(j,k);
                        for(ptrdiff_t i=j-1;i>=i0;i--)
                            x[i]-=x[j]*A[k-j+i];
                    }
                }
                else if(uplo=='L')
                {
                    for(size_t j=0;j<n;j++)
                    {
                        if(nounit)
                            x[j]=x[j]/A[0];
                        size_t im=minadd(j,k,n-1);
                        for(size_t i=j+1;i<=im;i++)
                            x[i]-=x[j]*A[i-j];
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
                        size_t i0=maxsub(j,k);
                        for(size_t i=i0;i<j;i++)
                            x[j]-=x[i]*A[k-j+i];
                        if(nounit)
                            x[j]=x[j]/A[k];
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        for(ptrdiff_t i=im;i>j;i--)
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
                        size_t i0=maxsub(j,k);
                        for(size_t i=i0;i<j;i++)
                            x[j]-=x[i]*conj(A[k-j+i]);
                        if(nounit)
                            x[j]=x[j]/conj(A[k]);
                        A+=ldA;
                    }
                }
                else if(uplo=='L')
                {
                    A+=n*ldA;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        for(ptrdiff_t i=im;i>j;i--)
                            x[j]-=x[i]*conj(A[i-j]);
                        if(nounit)
                            x[j]=x[j]/conj(A[0]);
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
                            x[jx]=x[jx]/A[k];
                        ptrdiff_t i0=maxsub(j,k);
                        for(ptrdiff_t i=j-1;i>=i0;i--)
                        {
                            ix-=incx;
                            x[ix]-=x[jx]*A[k-j+i];
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
                        size_t im=minadd(j,k,n-1);
                        for(size_t i=j+1;i<=im;i++)
                        {
                            ix+=incx;
                            x[ix]-=x[jx]*A[i-j];
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
                    size_t jx=kx;
                    for(size_t j=0;j<n;j++)
                    {
                        size_t i0=maxsub(j,k);
                        size_t ix=kx+i0*incx;
                        for(size_t i=i0;i<j;i++)
                        {
                            x[jx]-=x[ix]*A[k-j+i];
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/A[k];
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
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        ptrdiff_t ix=kx-(n-im)*incx;
                        for(ptrdiff_t i=im;i>j;i--)
                        {
                            x[jx]-=x[ix]*A[i-j];
                            ix-=incx;
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
                        size_t i0=maxsub(j,k);
                        size_t ix=kx+i0*incx;
                        for(size_t i=i0;i<j;i++)
                        {
                            x[jx]-=x[ix]*conj(A[k-j+i]);
                            ix+=incx;
                        }
                        if(nounit)
                            x[jx]=x[jx]/conj(A[k]);
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
                        A-=ldA;
                        ptrdiff_t im=minadd(j,k,n-1);
                        ptrdiff_t ix=kx-(n-im)*incx;
                        for(ptrdiff_t i=im;i>j;i--)
                        {
                            x[jx]-=x[ix]*conj(A[i-j]);
                            ix-=incx;
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
