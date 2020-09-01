#include "blas.h"
#include "trmv.h"
#include <cctype>

using std::toupper;
using tblas::trmv;

void strmv_(const char &Uplo, const char &Trans, const char &Diag, const int &n, float *A, const int &ldA, float *x, const int &incx)
{
    int info=0;
    char uplo=toupper(Uplo);
    char trans=toupper(Trans);
    char diag=toupper(Diag);
    if(trans=='C')
        trans='T';
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if((trans!='N')&&(trans!='T'))
        info=2;
    else if((diag!='U')&&(diag!='N'))
        info=3;
    else if(n<0)
        info=4;
    else if(ldA<n)
        info=6;
    else if(incx==0)
        info=8;
    if(info>0)
        xerbla_("STRMV ",info);
    else
        trmv(uplo,trans,diag,n,A,ldA,x,incx);
}
