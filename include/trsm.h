//
//  trsm.h
//
//  Purpose
//  =======
//
//  Solves one of the following systems:
//
//      A   * X = alpha * B  [trans='N' and side='L']
//
//      A^T * X = alpha * B  [trans='T' and side='L']
//
//      A^H * X = alpha * B  [trans='C' and side='L']
//
//      X * A   = alpha * B  [trans='N' and side='R']
//
//      X * A^T = alpha * B  [trans='T' and side='R']
//
//      X * A^H = alpha * B  [trans='C' and side='R']
//
//  where X and B are m-by-n matrices, and A is an n-by-n triangular matrix,
//  and alpha is a scalar.  If diag='U', the A is assumed to be unit triangular,
//  or diag='N' for nonunit triangular.
//
//  Arguments
//  ==========
//
//  side    specifies the side from which A is applied as above
//
//  uplo    specifies whether matrix is upper ('U') or lower ('L') triangular
//
//  trans   specifies the transpose operation for A as above
//
//  diag    specifies whether the matrix A is unit triangular ('U') or not ('N')
//
//  m       specifies the number of rows in the matrix B
//
//  n       specifies the order of the triangular matrix A
//
//  alpha   scalar multiple of the right-hand side B
//
//  A       triangular matrix stored in ldA-by-n array
//
//  ldA     column length of the matrix A, must be at least n
//
//  B       matrix of size m-by-n containing the right-hand side on entry,
//          and the solution X on exit, stored in ldB-by-b array
//
//  ldB     column length of the matrix B, must be at least m
//

#ifndef __trsm__
#define __trsm__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void trsm(char side, char uplo, char trans, char diag, size_t m, size_t n, T alpha, T *A, size_t ldA, T *B, size_t ldB)
    {
        const T zero(0.0);
        const T one(1.0);
        const bool nounit=(diag=='N');
        if((m==0)||(n==0))
            return;
        
        if(alpha==zero)
        {
            T *b=B;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    b[i]=zero;
                b+=ldB;
            }
        }
        else
        {
            if(side=='L')
            {
                if(trans=='N')
                {
                    if(uplo=='U')
                    {
                        T *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            T *a=A+m*ldA;
                            for(ptrdiff_t k=m-1;k>=0;k--)
                            {
                                a-=ldA;
                                if(nounit)
                                    b[k]=b[k]/a[k];
                                for(ptrdiff_t i=0;i<k;i++)
                                    b[i]-=b[k]*a[i];
                            }
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        T *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            T *a=A;
                            for(size_t k=0;k<m;k++)
                            {
                                if(nounit)
                                    b[k]=b[k]/a[k];
                                for(size_t i=k+1;i<m;i++)
                                    b[i]-=b[k]*a[i];
                                a+=ldA;
                            }
                            b+=ldB;
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        T *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            T *a=A;
                            for(size_t i=0;i<m;i++)
                            {
                                T t=alpha*b[i];
                                for(size_t k=0;k<i;k++)
                                    t-=a[k]*b[k];
                                if(nounit)
                                    t=t/a[i];
                                b[i]=t;
                                a+=ldA;
                            }
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        T *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            T *a=A+m*ldA;
                            for(ptrdiff_t i=m-1;i>=0;i--)
                            {
                                a-=ldA;
                                T t=alpha*b[i];
                                for(size_t k=i+1;k<m;k++)
                                    t-=a[k]*b[k];
                                if(nounit)
                                    t=t/a[i];
                                b[i]=t;
                            }
                            b+=ldB;
                        }
                    }
                }
            }
            else if(side=='R')
            {
                if(trans=='N')
                {
                    if(uplo=='U')
                    {
                        T *b=B;
                        T *a=A;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            T *bt=B;
                            for(size_t k=0;k<j;k++)
                            {
                                for(size_t i=0;i<m;i++)
                                    b[i]-=a[k]*bt[i];
                                bt+=ldB;
                            }
                            if(nounit)
                            {
                                T t=one/a[j];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            a+=ldA;
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        T *b=B+n*ldB;
                        T *a=A+n*ldA;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            a-=ldA;
                            b-=ldB;
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            T *bt=B+(j+1)*ldB;
                            for(size_t k=j+1;k<n;k++)
                            {
                                for(size_t i=0;i<m;i++)
                                    b[i]-=a[k]*bt[i];
                                bt+=ldB;
                            }
                            if(nounit)
                            {
                                T t=one/a[j];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        T *b=B+n*ldB;
                        T *a=A+n*ldA;
                        for(ptrdiff_t k=n-1;k>=0;k--)
                        {
                            a-=ldA;
                            b-=ldB;
                            if(nounit)
                            {
                                T t=one/a[k];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            T *bt=B;
                            for(ptrdiff_t j=0;j<k;j++)
                            {
                                T t=a[j];
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                        }
                    }
                    else if(uplo=='L')
                    {
                        T *a=A;
                        T *b=B;
                        for(size_t k=0;k<n;k++)
                        {
                            if(nounit)
                            {
                                T t=one/a[k];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            T *bt=B+(k+1)*ldB;
                            for(size_t j=k+1;j<n;j++)
                            {
                                T t=a[j];
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            a+=ldA;
                            b+=ldB;
                        }
                    }
                }
            }
        }
    }
    
    template <typename T>
    void trsm(char side, char uplo, char trans, char diag, size_t m, size_t n, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *B, size_t ldB)
    {
        const complex<T> zero(0.0);
        const complex<T> one(1.0);
        const bool nounit=(diag=='N');
        if((m==0)||(n==0))
            return;
        
        if(alpha==zero)
        {
            complex<T> *b=B;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    b[i]=zero;
                b+=ldB;
            }
        }
        else
        {
            if(side=='L')
            {
                if(trans=='N')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            complex<T> *a=A+m*ldA;
                            for(ptrdiff_t k=m-1;k>=0;k--)
                            {
                                a-=ldA;
                                if(nounit)
                                    b[k]=b[k]/a[k];
                                for(ptrdiff_t i=0;i<k;i++)
                                    b[i]-=b[k]*a[i];
                            }
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            complex<T> *a=A;
                            for(size_t k=0;k<m;k++)
                            {
                                if(nounit)
                                    b[k]=b[k]/a[k];
                                for(size_t i=k+1;i<m;i++)
                                    b[i]-=b[k]*a[i];
                                a+=ldA;
                            }
                            b+=ldB;
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            complex<T> *a=A;
                            for(size_t i=0;i<m;i++)
                            {
                                complex<T> t=alpha*b[i];
                                for(size_t k=0;k<i;k++)
                                    t-=a[k]*b[k];
                                if(nounit)
                                    t=t/a[i];
                                b[i]=t;
                                a+=ldA;
                            }
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            complex<T> *a=A+m*ldA;
                            for(ptrdiff_t i=m-1;i>=0;i--)
                            {
                                a-=ldA;
                                complex<T> t=alpha*b[i];
                                for(size_t k=i+1;k<m;k++)
                                    t-=a[k]*b[k];
                                if(nounit)
                                    t=t/a[i];
                                b[i]=t;
                            }
                            b+=ldB;
                        }
                    }
                }
                else if(trans=='C')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            complex<T> *a=A;
                            for(size_t i=0;i<m;i++)
                            {
                                complex<T> t=alpha*b[i];
                                for(size_t k=0;k<i;k++)
                                    t-=conj(a[k])*b[k];
                                if(nounit)
                                    t=t/conj(a[i]);
                                b[i]=t;
                                a+=ldA;
                            }
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<n;j++)
                        {
                            complex<T> *a=A+m*ldA;
                            for(ptrdiff_t i=m-1;i>=0;i--)
                            {
                                a-=ldA;
                                complex<T> t=alpha*b[i];
                                for(size_t k=i+1;k<m;k++)
                                    t-=conj(a[k])*b[k];
                                if(nounit)
                                    t=t/conj(a[i]);
                                b[i]=t;
                            }
                            b+=ldB;
                        }
                    }
                }
            }
            else if(side=='R')
            {
                if(trans=='N')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B;
                        complex<T> *a=A;
                        for(size_t j=0;j<n;j++)
                        {
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            complex<T> *bt=B;
                            for(size_t k=0;k<j;k++)
                            {
                                for(size_t i=0;i<m;i++)
                                    b[i]-=a[k]*bt[i];
                                bt+=ldB;
                            }
                            if(nounit)
                            {
                                complex<T> t=one/a[j];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            a+=ldA;
                            b+=ldB;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *b=B+n*ldB;
                        complex<T> *a=A+n*ldA;
                        for(ptrdiff_t j=n-1;j>=0;j--)
                        {
                            a-=ldA;
                            b-=ldB;
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            complex<T> *bt=B+(j+1)*ldB;
                            for(size_t k=j+1;k<n;k++)
                            {
                                for(size_t i=0;i<m;i++)
                                    b[i]-=a[k]*bt[i];
                                bt+=ldB;
                            }
                            if(nounit)
                            {
                                complex<T> t=one/a[j];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                        }
                    }
                }
                else if(trans=='T')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B+n*ldB;
                        complex<T> *a=A+n*ldA;
                        for(ptrdiff_t k=n-1;k>=0;k--)
                        {
                            a-=ldA;
                            b-=ldB;
                            if(nounit)
                            {
                                complex<T> t=one/a[k];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            complex<T> *bt=B;
                            for(ptrdiff_t j=0;j<k;j++)
                            {
                                complex<T> t=a[j];
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *a=A;
                        complex<T> *b=B;
                        for(size_t k=0;k<n;k++)
                        {
                            if(nounit)
                            {
                                complex<T> t=one/a[k];
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            complex<T> *bt=B+(k+1)*ldB;
                            for(size_t j=k+1;j<n;j++)
                            {
                                complex<T> t=a[j];
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            a+=ldA;
                            b+=ldB;
                        }
                    }
                }
                else if(trans=='C')
                {
                    if(uplo=='U')
                    {
                        complex<T> *b=B+n*ldB;
                        complex<T> *a=A+n*ldA;
                        for(ptrdiff_t k=n-1;k>=0;k--)
                        {
                            a-=ldA;
                            b-=ldB;
                            if(nounit)
                            {
                                complex<T> t=one/conj(a[k]);
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            complex<T> *bt=B;
                            for(ptrdiff_t j=0;j<k;j++)
                            {
                                complex<T> t=conj(a[j]);
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                        }
                    }
                    else if(uplo=='L')
                    {
                        complex<T> *a=A;
                        complex<T> *b=B;
                        for(size_t k=0;k<n;k++)
                        {
                            if(nounit)
                            {
                                complex<T> t=one/conj(a[k]);
                                for(size_t i=0;i<m;i++)
                                    b[i]*=t;
                            }
                            complex<T> *bt=B+(k+1)*ldB;
                            for(size_t j=k+1;j<n;j++)
                            {
                                complex<T> t=conj(a[j]);
                                for(size_t i=0;i<m;i++)
                                    bt[i]-=t*b[i];
                                bt+=ldB;
                            }
                            for(size_t i=0;i<m;i++)
                                b[i]*=alpha;
                            a+=ldA;
                            b+=ldB;
                        }
                    }
                }
            }
        }
    }
}
#endif
