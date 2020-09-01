#include "blas.h"
#include "symv.h"
#include <cctype>

using std::toupper;
using tblas::symv;

void ssymv_(const char &Uplo, const int &n, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy)
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
        xerbla_("SSYMV ",info);
    else
        symv(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
}
