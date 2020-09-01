#include "blas.h"
#include "hpr.h"
#include <cctype>

using std::toupper;
using tblas::hpr;

void zhpr_(const char &Uplo, const int &n, const double &alpha, complex<double> *x, const int &incx, complex<double> *A)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=5;
    if(info>0)
        xerbla_("ZHPR  ",info);
    else
        hpr(uplo,n,alpha,x,incx,A);
}
