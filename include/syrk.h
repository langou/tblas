//
//  syrk.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-matrix operations
//
//      C <- alpha * A * A^T + beta * C  [trans='N']
//
//      C <- alpha * A^T * A + beta * C  [trans='T']
//
//  where alpha and beta are scalars, C is a symmetric matrix,
//  and A is a n-by-k or k-by-n matrix.
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
//  alpha   scalar multiple of matrix-matrix product
//
//  A       matrix of size n-by-k if tran='N', or k-by-n if trans='T'
//
//  ldA     column length of the matrix A,
//          must be at least n if trans='N' or k if trans='T'
//
//  beta    scalar multiple of C
//
//  C       symmetric matrix of order n
//
//  ldC     column length of the matrix C, must be at least n
//

#ifndef __syrk__
#define __syrk__

#include <cstddef>

using std::size_t;

namespace tblas
{
    template <typename T>
    void syrk(char uplo, char trans, size_t n, size_t k, T alpha, T *A, size_t ldA, T beta, T *C, size_t ldC)
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
                    for(size_t l=0;l<k;l++)
                    {
                        T t=alpha*a[j];
                        for(size_t i=0;i<=j;i++)
                            c[i]+=t*a[i];
                        a+=ldA;
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
                    for(size_t l=0;l<k;l++)
                    {
                        T t=alpha*a[j];
                        for(size_t i=j;i<n;i++)
                            c[i]+=t*a[i];
                        a+=ldA;
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
                for(size_t j=0;j<n;j++)
                {
                    T *a=A;
                    for(size_t i=0;i<=j;i++)
                    {
                        T t=zero;
                        for(size_t l=0;l<k;l++)
                            t+=a[l]*at[l];
                        c[i]=alpha*t+beta*c[i];
                        a+=ldA;
                    }
                    at+=ldA;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                T *at=A;
                T *c=C;
                for(size_t j=0;j<n;j++)
                {
                    T *a=A+j*ldA;
                    for(size_t i=j;i<n;i++)
                    {
                        T t=zero;
                        for(size_t l=0;l<k;l++)
                            t+=a[l]*at[l];
                        c[i]=alpha*t+beta*c[i];
                        a+=ldA;
                    }
                    at+=ldA;
                    c+=ldC;
                }
            }
        }
    }
}
#endif
