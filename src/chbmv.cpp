#include "blas.h"
#include "hbmv.h"
#include <cctype>

using std::toupper;
using tblas::hbmv;

void chbmv_(const char &Uplo, const int &n, const int &k, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(k<0)
        info=3;
    else if(ldA<k+1)
        info=6;
    else if(incx==0)
        info=8;
    else if(incy==0)
        info=11;
    if(info>0)
        xerbla_("CHBMV ",info);
    else
        hbmv(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
}
