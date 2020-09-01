#include "blas.h"
#include "sbmv.h"
#include <cctype>

using std::toupper;
using tblas::sbmv;

void ssbmv_(const char &Uplo, const int &n, const int &k, const float &alpha, float *A, const int &ldA, float *x, const int &incx, const float &beta, float *y, const int &incy)
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
        xerbla_("SSBMV ",info);
    else
        sbmv(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
}
