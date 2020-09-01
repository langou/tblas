#include "blas.h"
#include "hpmv.h"
#include <cctype>

using std::toupper;
using tblas::hpmv;

void zhpmv_(const char &Uplo, const int &n, const complex<double> &alpha, complex<double> *A, complex<double> *x, const int &incx, const complex<double> &beta, complex<double> *y, const int &incy)
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
        xerbla_("ZHPMV ",info);
    else
        hpmv(uplo,n,alpha,A,x,incx,beta,y,incy);
}
