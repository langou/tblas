#include "blas.h"
#include "hemv.h"
#include <cctype>

using std::toupper;
using tblas::hemv;

void chemv_(const char &Uplo, const int &n, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *x, const int &incx, const complex<float> &beta, complex<float> *y, const int &incy)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(ldA<n)
        info=5;
    else if(incx==0)
        info=7;
    else if(incy==0)
        info=10;
    if(info>0)
        xerbla_("CHEMV ",info);
    else
        hemv(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
}
