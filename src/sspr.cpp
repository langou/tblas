#include "blas.h"
#include "spr.h"
#include <cctype>

using std::toupper;
using tblas::spr;

void sspr_(const char &Uplo, const int &n, const float &alpha, float *x, const int &incx, float *A)
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
        xerbla_("SSPR  ",info);
    else
        spr(uplo,n,alpha,x,incx,A);
}
