//
//  symm.h
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
//  where alpha and beta are scalars, A is a symmetric matrix,
//  and B and C are m-by-n matrices.
//
//  Arguments
//  =========
//
//  side    specifies the position of the symmetric matrix A above,
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

#ifndef __symm__
#define __symm__

#include <cstddef>

using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void symm(char side, char uplo, size_t m, size_t n, T alpha, T *A, size_t ldA, T *B, size_t ldB, T beta, T *C, size_t ldC)
    {
        const T zero(0.0);
        
        if(alpha==zero)
        {
            T *c=C;
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
                T *b=B;
                T *c=C;
                for(size_t j=0;j<n;j++)
                {
                    T *a=A;
                    for(size_t i=0;i<m;i++)
                    {
                        T t=alpha*b[i];
                        T s=zero;
                        for(size_t k=0;k<i;k++)
                        {
                            c[k]+=t*a[k];
                            s+=b[k]*a[k];
                        }
                        c[i]=beta*c[i]+t*a[i]+alpha*s;
                        a+=ldA;
                    }
                    b+=ldB;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                T *b=B;
                T *c=C;
                for(size_t j=0;j<n;j++)
                {
                    T *a=A+m*ldA;
                    for(ptrdiff_t i=m-1;i>=0;i--)
                    {
                        a-=ldA;
                        T t=alpha*b[i];
                        T s=zero;
                        for(size_t k=i+1;k<m;k++)
                        {
                            c[k]+=t*a[k];
                            s+=b[k]*a[k];
                        }
                        c[i]=beta*c[i]+t*a[i]+alpha*s;
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
                T *a=A;
                T *c=C;
                T *b=B;
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*a[j];
                    for(size_t i=0;i<m;i++)
                        c[i]=c[i]*beta+t*b[i];
                    T *at=A+(j+1)*ldA;
                    T *bt=B;
                    for(size_t k=0;k<j;k++)
                    {
                        T t=alpha*a[k];
                        for(size_t i=0;i<m;i++)
                            c[i]+=bt[i]*t;
                        bt+=ldB;
                    }
                    bt=B+(j+1)*ldB;
                    for(size_t k=j+1;k<n;k++)
                    {
                        T t=alpha*at[j];
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
                T *a=A;
                T *c=C;
                T *b=B;
                for(size_t j=0;j<n;j++)
                {
                    T t=alpha*a[j];
                    for(size_t i=0;i<m;i++)
                        c[i]=c[i]*beta+t*b[i];
                    T *bt=B;
                    T *at=A;
                    for(size_t k=0;k<j;k++)
                    {
                        T t=alpha*at[j];
                        for(size_t i=0;i<m;i++)
                            c[i]+=bt[i]*t;
                        at+=ldA;
                        bt+=ldB;
                    }
                    bt=B+(j+1)*ldB;
                    for(size_t k=j+1;k<n;k++)
                    {
                        T t=alpha*at[k];
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
