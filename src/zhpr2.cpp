#include "blas.h"
#include "hpr2.h"
#include <cctype>

using std::toupper;
using tblas::hpr2;

void zhpr2_(const char &Uplo, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=5;
    else if(incy==0)
        info=7;
    if(info>0)
        xerbla_("ZHPR2 ",info);
    else
        hpr2(uplo,n,alpha,x,incx,y,incy,A);
}
