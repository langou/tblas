//
//  herk.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-matrix operations
//
//      C <- alpha * A * A^H + beta * C  [trans='N']
//
//      C <- alpha * A^H * A + beta * C  [trans='C']
//
//  where alpha and beta are scalars, C is a Hermitian matrix,
//  and A is a complex n-by-k or k-by-n matrix.
//
//  Arguments
//  =========
//
//  uplo    specifies whether C is accessed as upper ('U') or lower ('L') triangular
//
//  trans   specifies whether the left matrices in the product are transposed,
//          must be either 'N' for no transpose, or 'C' for conjugate transpose
//
//  n       specifies the order of the Hermitian matrix C
//
//  k       specifies the inner dimension of the matrix-matrix products
//
//  alpha   scalar multiple of matrix-matrix product
//
//  A       complex matrix of size n-by-k if tran='N', or k-by-n if trans='C'
//
//  ldA     column length of the matrix A,
//          must be at least n if trans='N' or k if trans='C'
//
//  beta    scalar multiple of C
//
//  C       Hermitian matrix of order n
//
//  ldC     column length of the matrix C, must be at least n
//

#ifndef __herk__
#define __herk__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void herk(char uplo, char trans, size_t n, size_t k, T alpha, complex<T> *A, size_t ldA, T beta, complex<T> *C, size_t ldC)
    {
        const complex<T> zero(0.0);
        const T rzero(0.0);
        const T one(1.0);
        
        if((n==0)||(((alpha==rzero)||(k==0))&&(beta==one)))
            return;
        
        if(alpha==zero)
        {
            if(uplo=='U')
            {
                complex<T> *c=C;
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
                        for(size_t i=0;i<j;i++)
                            c[i]*=beta;
                        c[j]=beta*real(c[j]);
                        c+=ldC;
                    }
                }
            }
            else if(uplo=='L')
            {
                complex<T> *c=C;
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
                        c[j]=beta*real(c[j]);
                        for(size_t i=j+1;i<n;i++)
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
                complex<T> *c=C;
                for(size_t j=0;j<n;j++)
                {
                    for(size_t i=0;i<j;i++)
                        c[i]*=beta;
                    c[j]=beta*real(c[j]);
                    complex<T> *a=A;
                    for(size_t l=0;l<k;l++)
                    {
                        complex<T> t=alpha*conj(a[j]);
                        for(size_t i=0;i<j;i++)
                            c[i]+=t*a[i];
                        c[j]=real(c[j])+real(t*a[j]);
                        a+=ldA;
                    }
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                complex<T> *c=C;
                for(size_t j=0;j<n;j++)
                {
                    c[j]=beta*real(c[j]);
                    for(size_t i=j+1;i<n;i++)
                        c[i]*=beta;
                    complex<T> *a=A;
                    for(size_t l=0;l<k;l++)
                    {
                        complex<T> t=alpha*conj(a[j]);
                        c[j]=real(c[j])+real(t*a[j]);
                        for(size_t i=j+1;i<n;i++)
                            c[i]+=t*a[i];
                        a+=ldA;
                    }
                    c+=ldC;
                }
            }
        }
        else if(trans=='C')
        {
            if(uplo=='U')
            {
                complex<T> *c=C;
                complex<T> *at=A;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> *a=A;
                    for(size_t i=0;i<j;i++)
                    {
                        complex<T> t=zero;
                        for(size_t l=0;l<k;l++)
                            t+=conj(a[l])*at[l];
                        c[i]=alpha*t+beta*c[i];
                        a+=ldA;
                    }
                    T s=rzero;
                    for(size_t l=0;l<k;l++)
                        s+=real(conj(at[l])*at[l]);
                    c[j]=alpha*s+beta*real(c[j]);
                    at+=ldA;
                    c+=ldC;
                }
            }
            else if(uplo=='L')
            {
                complex<T> *at=A;
                complex<T> *c=C;
                for(size_t j=0;j<n;j++)
                {
                    T s=rzero;
                    for(size_t l=0;l<k;l++)
                        s+=real(conj(at[l])*at[l]);
                    c[j]=alpha*s+beta*real(c[j]);
                    complex<T> *a=A+(j+1)*ldA;
                    for(size_t i=j+1;i<n;i++)
                    {
                        complex<T> t=zero;
                        for(size_t l=0;l<k;l++)
                            t+=conj(a[l])*at[l];
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
