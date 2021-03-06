#include "blas.h"
#include "her2.h"
#include <cctype>

using std::toupper;
using tblas::her2;

void zher2_(const char &Uplo, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A, const int &ldA)
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
    else if(ldA<n)
        info=9;
    if(info>0)
        xerbla_("ZHER2 ",info);
    else
        her2(uplo,n,alpha,x,incx,y,incy,A,ldA);
}
