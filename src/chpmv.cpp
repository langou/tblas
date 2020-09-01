#include "blas.h"
#include "hpmv.h"
#include <cctype>

using std::toupper;
using tblas::hpmv;

void chpmv_(const char &Uplo, const int &n, const complex<float> &alpha, complex<float> *A, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=6;
    else if(incy==0)
        info=9;
    if(info>0)
        xerbla_("CHPMV ",info);
    else
        hpmv(uplo,n,alpha,A,x,incx,beta,y,incy);
}
