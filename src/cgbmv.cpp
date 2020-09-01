#include "blas.h"
#include "gbmv.h"
#include <cctype>

using std::toupper;
using tblas::gbmv;

void cgbmv_(const char &Trans, const int &m, const int &n, const int &kl, const int &ku, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy)
{
    int info=0;
    const char trans=toupper(Trans);
    if((trans!='N')&&(trans!='T')&&(trans!='C'))
        info=1;
    else if(m<0)
        info=2;
    else if(n<0)
        info=3;
    else if(kl<0)
        info=4;
    else if(ku<0)
        info=5;
    else if(ldA<kl+ku+1)
        info=8;
    else if(incx==0)
        info=10;
    else if(incy==0)
        info=13;
    if(info>0)
        xerbla_("CGBMV ",info);
    else
        gbmv(trans,m,n,kl,ku,alpha,A,ldA,x,incx,beta,y,incy);
}
