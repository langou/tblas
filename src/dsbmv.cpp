#include "blas.h"
#include "sbmv.h"
#include <cctype>

using std::toupper;
using tblas::sbmv;

void dsbmv_(const char &Uplo, const int &n, const int &k, const double &alpha, double *A, const int &ldA, double *x, const int &incx, const double &beta, double *y, const int &incy)
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
        xerbla_("DSBMV ",info);
    else
        sbmv(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
}
