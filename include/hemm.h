//
//  hemm.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-matrix operations
//
//      C <- alpha * A * B + beta * C  [side='L']
//
//      C <- alpha * B * A + beta * C  [side='R']
//
//  where alpha and beta are scalars, A is a Hermitian matrix,
//  and B and C are m-by-n matrices.
//
//  Arguments
//  =========
//
//  side    specifies the position of the Hermitian matrix A above,
//          must be either 'L' or 'R'
//
//  uplo    specifies whether A is accessed as upper ('U') or lower ('L') triangular
//
//  m       number of rows of the matrices B and C
//
//  n       number of columns of the matrices B and C
//
//  alpha   scalar multiple of A
//
//  A       triangular matrix of order m if side='L' or order n if side='R'
//
//  ldA     column length of the matrix A, must be at least m if side='L',
//          or at least n if side='R'
//
//  B       matrix of size m-by-n
//
//  ldB     column length of the matrix B, must be at least m
//
//  beta    scalar multiple of C
//
//  C       matrix of size m-by-n
//
//  ldC     column length of the matrix C, must be at least m
//

#ifndef __hemm__
#define __hemm__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void hemm(char side, char uplo, size_t m, size_t n, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *B, size_t ldB, complex<T> beta, complex<T> *C, size_t ldC)
    {
        const complex<T> zero(0.0);
        const complex<T> one(1.0);
        
        if((m==0)||(n==0)||((alpha==zero)&&(beta==one)))
            return;
        
        if(alpha==zero)
        {
            complex<T> *c=C;
            if(beta==zero)
            {
                for(size_t j=0;j<n;j++)
                {
                    for(size_t i=0;i<m;i++)
                        c[i]=zero;
                    c+=ldC;
                }
            }
            else
            {
                for(size_t j=0;j<n;j++)
                {
                    for(size_t i=0;i<m;i++)
                        c[i]*=beta;
                    c+=ldC;
                }
            }
        }
        else if(side=='L')
        {
            if(uplo=='U')
            {
                complex<T> *b=B;
                complex<T> *c=C;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> *a=A;
                    for(size_t i=0;i<m;i++)
                    {
                        complex<T> t=alpha*b[i];
                        complex<T> s=zero;
                        for(size_t k=0;k<i;k++)
                        {
                            c[k]+=t*a[k];
                            s+=b[k]*conj(a[k]);
                        }
                        c[i]=beta*c[i]+t*real(a[i])+alpha*s;
                        a+=ldA;
                    }
                    b+=ldB;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                complex<T> *b=B;
                complex<T> *c=C;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> *a=A+m*ldA;
                    for(ptrdiff_t i=m-1;i>=0;i--)
                    {
                        a-=ldA;
                        complex<T> t=alpha*b[i];
                        complex<T> s=zero;
                        for(size_t k=i+1;k<m;k++)
                        {
                            c[k]+=t*a[k];
                            s+=b[k]*conj(a[k]);
                        }
                        c[i]=beta*c[i]+t*real(a[i])+alpha*s;
                    }
                    b+=ldB;
                    c+=ldC;
                }
            }
        }
        else if(side=='R')
        {
            if(uplo=='U')
            {
                complex<T> *a=A;
                complex<T> *c=C;
                complex<T> *b=B;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*real(a[j]);
                    for(size_t i=0;i<m;i++)
                        c[i]=c[i]*beta+t*b[i];
                    complex<T> *at=A+(j+1)*ldA;
                    complex<T> *bt=B;
                    for(size_t k=0;k<j;k++)
                    {
                        t=alpha*a[k];
                        for(size_t i=0;i<m;i++)
                            c[i]+=bt[i]*t;
                        bt+=ldB;
                    }
                    bt=B+(j+1)*ldB;
                    for(size_t k=j+1;k<n;k++)
                    {
                        complex<T> t=alpha*conj(at[j]);
                        for(size_t i=0;i<m;i++)
                            c[i]+=t*bt[i];
                        at+=ldA;
                        bt+=ldB;
                    }
                    a+=ldA;
                    b+=ldB;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                complex<T> *a=A;
                complex<T> *c=C;
                complex<T> *b=B;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> t=alpha*real(a[j]);
                    for(size_t i=0;i<m;i++)
                        c[i]=c[i]*beta+t*b[i];
                    complex<T> *bt=B;
                    complex<T> *at=A;
                    for(size_t k=0;k<j;k++)
                    {
                        complex<T> t=alpha*conj(at[j]);
                        for(size_t i=0;i<m;i++)
                            c[i]+=bt[i]*t;
                        at+=ldA;
                        bt+=ldB;
                    }
                    bt=B+(j+1)*ldB;
                    for(size_t k=j+1;k<n;k++)
                    {
                        complex<T> t=alpha*at[k];
                        for(size_t i=0;i<m;i++)
                            c[i]+=t*bt[i];
                        bt+=ldB;
                    }
                    a+=ldA;
                    b+=ldB;
                    c+=ldC;
                }
            }
        }
    }
}
#endif
