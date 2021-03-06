#include "blas.h"
#include "gemv.h"
#include <cctype>

using std::toupper;
using tblas::gemv;

void cgemv_(const char &Trans, const int &m, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy)
{
    int info=0;
    char trans=toupper(Trans);
    if((trans!='N')&&(trans!='T')&&(trans!='C'))
        info=1;
    else if(m<0)
        info=2;
    else if(n<0)
        info=3;
    else if(ldA<m)
        info=6;
    else if(incx==0)
        info=8;
    else if(incy==0)
        info=11;
    if(info>0)
        xerbla_("CGEMV ",info);
    else
        gemv(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy);
}
