//
//  gemm.h
//
//  Purpose
//  =======
//
//  Performs one of the following matrix-matrix operations
//
//      C <- alpha * A   * B   + beta * C  [transA='N' and transB='N']
//
//      C <- alpha * A^T * B   + beta * C  [transA='T' and transB='N']
//
//      C <- alpha * A   * B^T + beta * C  [transA='N' and transB='T']
//
//      C <- alpha * A^T * B^T + beta * C  [transA='T' and transB='T']
//
//      C <- alpha * A^H * B   + beta * C  [transA='C' and transB='N']
//
//      C <- alpha * A   * B^H + beta * C  [transA='N' and transB='C']
//
//      C <- alpha * A^H * B^H + beta * C  [transA='C' and transB='C']
//
//      C <- alpha * A^T * B^H + beta * C  [transA='T' and transB='C']
//
//      C <- alpha * A^H * B^T + beta * C  [transA='C' and transB='T']
//
//  where alpha and beta are scalars, and A, B and C are matrices where
//  C is m-by-n.
//
//  Arguments
//  =========
//
//  transA  specifies the transpose operation for A as above
//
//  transB  specifies the transpose operation for B as above
//
//  m       number of rows of the matrix C
//
//  n       number of columns of the matrix C
//
//  k       inner matrix product dimension
//
//  alpha   scalar multiple of A
//
//  A       matrix of size m-by-k if transA='N', or k-by-m otherwise
//
//  ldA     column length of the matrix A, must be at least m if transA='N',
//          or at least k otherwise
//
//  B       matrix of size k-by-n if transB='N', or n-by-k otherwise
//
//  ldB     column length of the matrix B, must be ast least k if transB='N',
//          or at least n otherwise
//
//  beta    scalar multiple of C
//
//  C       matrix of size m-by-n
//
//  ldC     column length of the matrix C, must be at least m
//

#ifndef __gemm__
#define __gemm__

#include <cstddef>
#include <complex>

using std::size_t;
using std::complex;

namespace tblas
{
    template <typename T>
    void gemm(char transA, char transB, size_t m, size_t n, size_t k, T alpha, T *A, size_t ldA, T *B, size_t ldB, T beta, T *C, size_t ldC)
    {
        const T zero(0.0);
        const T one(1.0);
        
        if((m==0)||(n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
            return;
        
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
        else if((transA=='N')&&(transB=='N'))
        {
            T *c=C;
            T *b=B;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    c[i]*=beta;
                T *a=A;
                for(size_t l=0;l<k;l++)
                {
                    for(size_t i=0;i<m;i++)
                    {
                        c[i]+=alpha*b[l]*a[i];
                    }
                    a+=ldA;
                }
                b+=ldB;
                c+=ldC;
            }
        }
        else if((transA=='T')&&(transB=='N'))
        {
            T *b=B;
            T *c=C;
            for(size_t j=0;j<n;j++)
            {
                T *a=A;
                for(size_t i=0;i<m;i++)
                {
                    T s=zero;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=a[l]*b[l];
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                b+=ldB;
                c+=ldC;
            }
        }
        else if((transA=='N')&&(transB=='T'))
        {
            T *c=C;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    c[i]*=beta;
                T *a=A;
                T *b=B;
                for(size_t l=0;l<k;l++)
                {
                    T t=alpha*b[j];
                    for(size_t i=0;i<m;i++)
                    {
                        c[i]+=t*a[i];
                    }
                    a+=ldA;
                    b+=ldB;
                }
                c+=ldC;
            }
        }
        else if((transA=='T')&&(transB=='T'))
        {
            T *c=C;
            for(size_t j=0;j<n;j++)
            {
                T *a=A;
                for(size_t i=0;i<m;i++)
                {
                    T s=zero;
                    T *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=a[l]*b[j];
                        b+=ldB;
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                c+=ldC;
            }
            
        }
    }

    template <typename T>
    void gemm(char transA, char transB, size_t m, size_t n, size_t k, complex<T> alpha, complex<T> *A, size_t ldA, complex<T> *B, size_t ldB, complex<T> beta, complex<T> *C, size_t ldC)
    {
        const complex<T> zero(0.0);
        const complex<T> one(1.0);
        
        if((m==0)||(n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
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
        else if((transA=='N')&&(transB=='N'))
        {
            complex<T> *c=C;
            complex<T> *b=B;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    c[i]*=beta;
                complex<T> *a=A;
                for(size_t l=0;l<k;l++)
                {
                    for(size_t i=0;i<m;i++)
                    {
                        c[i]+=alpha*b[l]*a[i];
                    }
                    a+=ldA;
                }
                b+=ldB;
                c+=ldC;
            }
        }
        else if((transA=='T')&&(transB=='N'))
        {
            complex<T> *b=B;
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=a[l]*b[l];
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                b+=ldB;
                c+=ldC;
            }
        }
        else if((transA=='C')&&(transB=='N'))
        {
            complex<T> *b=B;
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=conj(a[l])*b[l];
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                b+=ldB;
                c+=ldC;
            }
        }
        else if((transA=='N')&&(transB=='T'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    c[i]*=beta;
                complex<T> *a=A;
                complex<T> *b=B;
                for(size_t l=0;l<k;l++)
                {
                    complex<T> t=alpha*b[j];
                    for(size_t i=0;i<m;i++)
                    {
                        c[i]+=t*a[i];
                    }
                    a+=ldA;
                    b+=ldB;
                }
                c+=ldC;
            }
        }
        else if((transA=='N')&&(transB=='C'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                for(size_t i=0;i<m;i++)
                    c[i]*=beta;
                complex<T> *a=A;
                complex<T> *b=B;
                for(size_t l=0;l<k;l++)
                {
                    complex<T> t=alpha*conj(b[j]);
                    for(size_t i=0;i<m;i++)
                    {
                        c[i]+=t*a[i];
                    }
                    a+=ldA;
                    b+=ldB;
                }
                c+=ldC;
            }
        }
        else if((transA=='T')&&(transB=='T'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=a[l]*b[j];
                        b+=ldB;
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                c+=ldC;
            }
        }
        else if((transA=='C')&&(transB=='T'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=conj(a[l])*b[j];
                        b+=ldB;
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                c+=ldC;
            }
        }
        else if((transA=='T')&&(transB=='C'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=a[l]*conj(b[j]);
                        b+=ldB;
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                c+=ldC;
            }
        }
        else if((transA=='C')&&(transB=='C'))
        {
            complex<T> *c=C;
            for(size_t j=0;j<n;j++)
            {
                complex<T> *a=A;
                for(size_t i=0;i<m;i++)
                {
                    complex<T> s=zero;
                    complex<T> *b=B;
                    for(size_t l=0;l<k;l++)
                    {
                        s+=conj(a[l])*conj(b[j]);
                        b+=ldB;
                    }
                    c[i]=alpha*s+beta*c[i];
                    a+=ldA;
                }
                c+=ldC;
            }
        }
    }
}
#endif
