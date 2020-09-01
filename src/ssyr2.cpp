#include "blas.h"
#include "syr2.h"
#include <cctype>

using std::toupper;
using tblas::syr2;

void ssyr2_(const char &Uplo, const int &n, const float &alpha, float *x, const int &incx, float *y, const int &incy, float *A, const int &ldA)
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
        xerbla_("SSYR2 ",info);
    else
        syr2(uplo,n,alpha,x,incx,y,incy,A,ldA);
}
