#include "blas.h"
#include "tbmv.h"
#include <cctype>

using std::toupper;
using tblas::tbmv;

void dtbmv_(const char &Uplo, const char &Trans, const char &Diag, const int &n, const int &k, double *A, const int &ldA, double *x, const int &incx)
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
    else if(k<0)
        info=5;
    else if(ldA<k+1)
        info=7;
    else if(incx==0)
        info=9;
    if(info>0)
        xerbla_("DTBMV ",info);
    else
        tbmv(uplo,trans,diag,n,k,A,ldA,x,incx);
}
