//
//  syr2k.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-matrix operations
//
//      C <- alpha * A * B^H + alpha * B * A^T + beta * C [trans='N']
//
//      C <- alpha * A^T * B + alpha * B^T * A + beta * C [trans='T']
//
//  where alpha and beta are scalars, C is a symmetric matrix,
//  and A and B are n-by-k or k-by-n matrices.
//
//  Arguments
//  =========
//
//  uplo    specifies whether C is accessed as upper ('U') or lower ('L') triangular
//
//  trans   specifies whether the left matrices in the product are transposed,
//          must be either 'N' for no transpose, or 'T' for transpose
//
//  n       specifies the order of the symmetric matrix C
//
//  k       specifies the inner dimension of the matrix-matrix products
//
//  alpha   scalar multiple of the matrix-matrix product
//
//  A       matrix of size n-by-k if tran='N', or k-by-n if trans='T'
//
//  ldA     column length of the matrix A, must be at least n if trans='N' or k if trans='T'
//
//  B       matrix of size n-by-k if tran='N', or k-by-n if trans='T'
//
//  ldB     column length of the matrix B, must be at least n if trans='N' or k if trans='T'
//
//  beta    scalar multiple of C
//
//  C       symmetric matrix of order n
//
//  ldC     column length of the matrix C, must be at least n
//

#ifndef __syr2k__
#define __syr2k__

#include <cstddef>

using std::size_t;

namespace tblas
{
    template <typename T>
    void syr2k(char uplo, char trans, size_t n, size_t k, T alpha, T *A, size_t ldA, T *B, size_t ldB, T beta, T *C, size_t ldC)
    {
        const T zero(0.0);
        const T one(1.0);
        
        if((n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
            return;
        
        if(alpha==zero)
        {
            if(uplo=='U')
            {
                T *c=C;
                if(beta==zero)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<=j;i++)
                            c[i]=zero;
                        c+=ldC;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=0;i<=j;i++)
                            c[i]*=beta;
                        c+=ldC;
                    }
                }
            }
            else if(uplo=='L')
            {
                T *c=C;
                if(beta==zero)
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=j;i<n;i++)
                            c[i]=zero;
                        c+=ldC;
                    }
                }
                else
                {
                    for(size_t j=0;j<n;j++)
                    {
                        for(size_t i=j;i<n;i++)
                            c[i]*=beta;
                        c+=ldC;
                    }
                }
            }
        }
        else if(trans=='N')
        {
            if(uplo=='U')
            {
                T *c=C;
                for(size_t j=0;j<n;j++)
                {
                    for(size_t i=0;i<=j;i++)
                        c[i]*=beta;
                    T *a=A;
                    T *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        T s=alpha*b[j];
                        T t=alpha*a[j];
                        for(size_t i=0;i<=j;i++)
                            c[i]+=s*a[i]+t*b[i];
                        a+=ldA;
                        b+=ldB;
                    }
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                T *c=C;
                for(size_t j=0;j<n;j++)
                {
                    for(size_t i=j;i<n;i++)
                        c[i]*=beta;
                    T *a=A;
                    T *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        T s=alpha*b[j];
                        T t=alpha*a[j];
                        for(size_t i=j;i<n;i++)
                            c[i]+=s*a[i]+t*b[i];
                        a+=ldA;
                        b+=ldB;
                    }
                    c+=ldC;
                }
            }
        }
        else if(trans=='T')
        {
            if(uplo=='U')
            {
                T *c=C;
                T *at=A;
                T *bt=B;
                for(size_t j=0;j<n;j++)
                {
                    T *a=A;
                    T *b=B;
                    for(size_t i=0;i<=j;i++)
                    {
                        T s=zero;
                        T t=zero;
                        for(size_t l=0;l<k;l++)
                        {
                            s+=b[l]*at[l];
                            t+=a[l]*bt[l];
                        }
                        c[i]=alpha*t+alpha*s+beta*c[i];
                        a+=ldA;
                        b+=ldB;
                    }
                    at+=ldA;
                    bt+=ldB;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                T *c=C;
                T *at=A;
                T *bt=B;
                for(size_t j=0;j<n;j++)
                {
                    T *a=A+j*ldA;
                    T *b=B+j*ldB;
                    for(size_t i=j;i<n;i++)
                    {
                        T s=zero;
                        T t=zero;
                        for(size_t l=0;l<k;l++)
                        {
                            s+=b[l]*at[l];
                            t+=a[l]*bt[l];
                        }
                        c[i]=alpha*t+alpha*s+beta*c[i];
                        a+=ldA;
                        b+=ldB;
                    }
                    at+=ldA;
                    bt+=ldB;
                    c+=ldC;
                }
            }
        }
    }
}
#endif
