#include "blas.h"
#include "hemv.h"
#include <cctype>

using std::toupper;
using tblas::hemv;

void zhemv_(const char &Uplo, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy)
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
        xerbla_("ZHEMV ",info);
    else
        hemv(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
}
