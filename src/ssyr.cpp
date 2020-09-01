#include "blas.h"
#include "syr.h"
#include <cctype>

using std::toupper;
using tblas::syr;

void ssyr_(const char &Uplo, const int &n, const float &alpha, float *x, const int &incx, float *A, const int &ldA)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=5;
    else if(ldA<n)
        info=7;
    if(info>0)
        xerbla_("SSYR  ",info);
    else
        syr(uplo,n,alpha,x,incx,A,ldA);
}
