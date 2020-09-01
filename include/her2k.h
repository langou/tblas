//
//  her2k.h
//
//  Purpose
//  =======
//
//  Performs one of the matrix-matrix operations
//
//      C <- alpha * A * B^H + conj(alpha) * B * A^H + beta * C [trans='N']
//
//      C <- alpha * A^H * B + conj(alpha) * B^H * A + beta * C [trans='C']
//
//  where alpha and beta are scalars, C is a Hermitian matrix,
//  and A and B are complex n-by-k or k-by-n matrices.
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
//  alpha   scalar multiple of the matrix-matrix product
//
//  A       complex matrix of size n-by-k if tran='N', or k-by-n if trans='C'
//
//  ldA     column length of the matrix A,
//          must be at least n if trans='N' or k if trans='C'
//
//  B       complex matrix of size n-by-k if tran='N', or k-by-n if trans='C'
//
//  ldB     column length of the matrix B,
//          must be at least n if trans='N' or k if trans='C'
//
//  beta    scalar multiple of C
//
//  C       Hermitian matrix of order n
//
//  ldC     column length of the matrix C, must be at least n
//

#ifndef __her2k__
#define __her2k__

#include <complex>
#include <cstddef>

using std::complex;
using std::size_t;
using std::ptrdiff_t;

namespace tblas
{
    template <typename T>
    void her2k(char uplo, char trans, size_t n, size_t k, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *B, size_t ldB, T beta, complex<T> *C, size_t ldC)
    {
        const complex<T> zero(0.0);
        const complex<T> one(1.0);
        
        if((n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
            return;
        
        if(alpha==zero)
        {
            complex<T> *c=C;
            if(uplo=='U')
            {
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
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        complex<T> s=alpha*conj(b[j]);
                        complex<T> t=conj(alpha*a[j]);
                        for(size_t i=0;i<j;i++)
                            c[i]+=s*a[i]+t*b[i];
                        c[j]=real(c[j])+real(a[j]*s+b[j]*t);
                        a+=ldA;
                        b+=ldB;
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
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        complex<T> s=alpha*conj(b[j]);
                        complex<T> t=conj(alpha*a[j]);
                        c[j]=real(c[j])+real(s*a[j]+t*b[j]);
                        for(size_t i=j+1;i<n;i++)
                            c[i]+=s*a[i]+t*b[i];
                        a+=ldA;
                        b+=ldB;
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
                complex<T> *bt=B;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> *a=A;
                    complex<T> *b=B;
                    for(size_t i=0;i<=j;i++)
                    {
                        complex<T> s=zero;
                        complex<T> t=zero;
                        for(size_t l=0;l<k;l++)
                        {
                            s+=conj(b[l])*at[l];
                            t+=conj(a[l])*bt[l];
                        }
                        if(i<j)
                            c[i]=alpha*t+conj(alpha)*s+beta*c[i];
                        else
                            c[j]=real(alpha*t+conj(alpha)*s)+beta*real(c[j]);
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
                complex<T> *c=C;
                complex<T> *at=A;
                complex<T> *bt=B;
                for(size_t j=0;j<n;j++)
                {
                    complex<T> *a=A+j*ldA;
                    complex<T> *b=B+j*ldB;
                    for(size_t i=j;i<n;i++)
                    {
                        complex<T> s=zero;
                        complex<T> t=zero;
                        for(size_t l=0;l<k;l++)
                        {
                            s+=conj(b[l])*at[l];
                            t+=conj(a[l])*bt[l];
                        }
                        if(i>j)
                            c[i]=alpha*t+conj(alpha)*s+beta*c[i];
                        else
                            c[j]=real(alpha*t+conj(alpha)*s)+beta*real(c[j]);
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
