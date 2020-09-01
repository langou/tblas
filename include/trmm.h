//
//  trmm.h
//
//  Purpose
//  =======
//
//  Performs multiplications of a triangular matrix with a rectangular matrix:
//
//      B <- alpha * A   * B  [trans='N' and side='L']
//
//      B <- alpha * A^T * B  [trans='T' and side='L']
//
//      B <- alpha * A^H * B  [trans='C' and side='L']
//
//      B <- alpha * B * A    [trans='N' and side='R']
//
//      B <- alpha * B * A^T  [trans='T' and side='R']
//
//      B <- alpha * B * A^H  [trans='C' and side='R']
//
//  where B is a m-by-n matrix, A is n-by-n triangular matrix and alpha is a scalar.
//  If diag='U' A is assumed to be unit triangular, or diag='N' for nonunit triangular.
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
//  B       matrix of size m-by-n stored in ldB-by-b array
//
//  ldB     column length of the matrix B, must be at least m
//

#ifndef __trmm__
#define __trmm__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void trmm(char side, char uplo, char trans, char diag, size_t m, size_t n, T alpha, T *A, size_t ldA, T *B, size_t ldB)
    {
        const T zero(0.0);
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
        else if(side=='L')
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    T *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        T *a=A;
                        for(size_t k=0;k<m;k++)
                        {
                            T t=alpha*b[k];
                            for(size_t i=0;i<k;i++)
                                b[i]+=t*a[i];
                            if(nounit)
                                t*=a[k];
                            b[k]=t;
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
                        for(ptrdiff_t k=m-1;k>=0;k--)
                        {
                            a-=ldA;
                            T t=alpha*b[k];
                            b[k]=t;
                            if(nounit)
                                b[k]*=a[k];
                            for(size_t i=k+1;i<m;i++)
                                b[i]+=t*a[i];
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
                        T *a=A+m*ldA;
                        for(ptrdiff_t i=m-1;i>=0;i--)
                        {
                            a-=ldA;
                            T t=b[i];
                            if(nounit)
                                t*=a[i];
                            for(ptrdiff_t k=0;k<i;k++)
                                t+=a[k]*b[k];
                            b[i]=alpha*t;
                        }
                        b+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    T *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        T *a=A;
                        for(size_t i=0;i<m;i++)
                        {
                            T t=b[i];
                            if(nounit)
                                t*=a[i];
                            for(size_t k=i+1;k<m;k++)
                                t+=a[k]*b[k];
                            b[i]=alpha*t;
                            a+=ldA;
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
                    T *a=A+n*ldA;
                    T *b=B+n*ldB;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        a-=ldA;
                        b-=ldB;
                        T s=(nounit)?alpha*a[j]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=s;
                        T *bt=B;
                        for(ptrdiff_t k=0;k<j;k++)
                        {
                            T t=alpha*a[k];
                            for(size_t i=0;i<m;i++)
                                b[i]+=t*bt[i];
                            bt+=ldB;
                        }
                    }
                }
                else if(uplo=='L')
                {
                    T *a=A;
                    T *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        T s=(nounit)?alpha*a[j]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=s;
                        T *bt=b+ldB;
                        for(size_t k=j+1;k<n;k++)
                        {
                            T t=alpha*a[k];
                            for(size_t i=0;i<m;i++)
                                b[i]+=t*bt[i];
                            bt+=ldB;
                        }
                        a+=ldA;
                        b+=ldB;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    T *a=A;
                    T *bt=B;
                    for(size_t k=0;k<n;k++)
                    {
                        T *b=B;
                        for(size_t j=0;j<k;j++)
                        {
                            T t=alpha*a[j];
                            for(size_t i=0;i<m;i++)
                                b[i]+=t*bt[i];
                            b+=ldB;
                        }
                        T s=(nounit)?alpha*a[k]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=s;
                        a+=ldA;
                        bt+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    T *a=A+n*ldA;
                    T *bt=B+n*ldB;
                    for(ptrdiff_t k=n-1;k>=0;k--)
                    {
                        a-=ldA;
                        bt-=ldB;
                        T *b=B+(k+1)*ldB;
                        for(size_t j=k+1;j<n;j++)
                        {
                            T t=alpha*a[j];
                            for(size_t i=0;i<m;i++)
                                b[i]+=t*bt[i];
                            b+=ldB;
                        }
                        T s=(nounit)?alpha*a[k]:alpha;
                        for(size_t i=0;i<m;i++)
                            bt[i]*=s;
                    }
                }
            }
        }
    }

    template <typename T>
    void trmm(char side, char uplo, char trans, char diag, size_t m, size_t n, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *B, size_t ldB)
    {
        const complex<T> zero(0.0);
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
        else if(side=='L')
        {
            if(trans=='N')
            {
                if(uplo=='U')
                {
                    complex<T> *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> *a=A;
                        for(size_t k=0;k<m;k++)
                        {
                            complex<T> t=alpha*b[k];
                            for(size_t i=0;i<k;i++)
                                b[i]+=t*a[i];
                            if(nounit)
                                t*=a[k];
                            b[k]=t;
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
                        for(ptrdiff_t k=m-1;k>=0;k--)
                        {
                            a-=ldA;
                            complex<T> t=alpha*b[k];
                            b[k]=t;
                            if(nounit)
                                b[k]*=a[k];
                            for(size_t i=k+1;i<m;i++)
                                b[i]+=t*a[i];
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
                        complex<T> *a=A+m*ldA;
                        for(ptrdiff_t i=m-1;i>=0;i--)
                        {
                            a-=ldA;
                            complex<T> t=b[i];
                            if(nounit)
                                t*=a[i];
                            for(ptrdiff_t k=0;k<i;k++)
                                t+=a[k]*b[k];
                            b[i]=alpha*t;
                        }
                        b+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    complex<T> *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> *a=A;
                        for(size_t i=0;i<m;i++)
                        {
                            complex<T> t=b[i];
                            if(nounit)
                                t*=a[i];
                            for(size_t k=i+1;k<m;k++)
                                t+=a[k]*b[k];
                            b[i]=alpha*t;
                            a+=ldA;
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
                        complex<T> *a=A+m*ldA;
                        for(ptrdiff_t i=m-1;i>=0;i--)
                        {
                            a-=ldA;
                            complex<T> t=b[i];
                            if(nounit)
                                t*=conj(a[i]);
                            for(ptrdiff_t k=0;k<i;k++)
                                t+=conj(a[k])*b[k];
                            b[i]=alpha*t;
                        }
                        b+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    complex<T> *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> *a=A;
                        for(size_t i=0;i<m;i++)
                        {
                            complex<T> t=b[i];
                            if(nounit)
                                t*=conj(a[i]);
                            for(size_t k=i+1;k<m;k++)
                                t+=conj(a[k])*b[k];
                            b[i]=alpha*t;
                            a+=ldA;
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
                    complex<T> *a=A+n*ldA;
                    complex<T> *b=B+n*ldB;
                    for(ptrdiff_t j=n-1;j>=0;j--)
                    {
                        a-=ldA;
                        b-=ldB;
                        complex<T> t=(nounit)?alpha*a[j]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=t;
                        complex<T> *bt=B;
                        for(ptrdiff_t k=0;k<j;k++)
                        {
                            complex<T> s=alpha*a[k];
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            bt+=ldB;
                        }
                    }
                }
                else if(uplo=='L')
                {
                    complex<T> *a=A;
                    complex<T> *b=B;
                    for(size_t j=0;j<n;j++)
                    {
                        complex<T> t=(nounit)?alpha*a[j]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=t;
                        complex<T> *bt=b+ldB;
                        for(size_t k=j+1;k<n;k++)
                        {
                            complex<T> s=alpha*a[k];
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            bt+=ldB;
                        }
                        a+=ldA;
                        b+=ldB;
                    }
                }
            }
            else if(trans=='T')
            {
                if(uplo=='U')
                {
                    complex<T> *a=A;
                    complex<T> *bt=B;
                    for(size_t k=0;k<n;k++)
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<k;j++)
                        {
                            complex<T> s=alpha*a[j];
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            b+=ldB;
                        }
                        complex<T> t=(nounit)?alpha*a[k]:alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=t;
                        a+=ldA;
                        bt+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    complex<T> *a=A+n*ldA;
                    complex<T> *bt=B+n*ldB;
                    for(ptrdiff_t k=n-1;k>=0;k--)
                    {
                        a-=ldA;
                        bt-=ldB;
                        complex<T> *b=B+(k+1)*ldB;
                        for(size_t j=k+1;j<n;j++)
                        {
                            complex<T> s=alpha*a[j];
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            b+=ldB;
                        }
                        complex<T> t=(nounit)?alpha*a[k]:alpha;
                        for(size_t i=0;i<m;i++)
                            bt[i]*=t;
                    }
                }
            }
            else if(trans=='C')
            {
                if(uplo=='U')
                {
                    complex<T> *a=A;
                    complex<T> *bt=B;
                    for(size_t k=0;k<n;k++)
                    {
                        complex<T> *b=B;
                        for(size_t j=0;j<k;j++)
                        {
                            complex<T> s=alpha*conj(a[j]);
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            b+=ldB;
                        }
                        complex<T> t=(nounit)?alpha*conj(a[k]):alpha;
                        for(size_t i=0;i<m;i++)
                            b[i]*=t;
                        a+=ldA;
                        bt+=ldB;
                    }
                }
                else if(uplo=='L')
                {
                    complex<T> *a=A+n*ldA;
                    complex<T> *bt=B+n*ldB;
                    for(ptrdiff_t k=n-1;k>=0;k--)
                    {
                        a-=ldA;
                        bt-=ldB;
                        complex<T> *b=B+(k+1)*ldB;
                        for(size_t j=k+1;j<n;j++)
                        {
                            complex<T> s=alpha*conj(a[j]);
                            for(size_t i=0;i<m;i++)
                                b[i]+=s*bt[i];
                            b+=ldB;
                        }
                        complex<T> t=(nounit)?alpha*conj(a[k]):alpha;
                        for(size_t i=0;i<m;i++)
                            bt[i]*=t;
                    }
                }
            }
        }
    }
    
}
#endif
