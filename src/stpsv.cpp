#include "blas.h"
#include "tpsv.h"
#include <cctype>

using std::toupper;
using tblas::tpsv;

void stpsv_(const char &Uplo, const char &Trans, const char &Diag, const int &n, float *A, float *x, const int &incx)
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
    else if(incx==0)
        info=7;
    if(info>0)
        xerbla_("STPSV ",info);
    else
        tpsv(uplo,trans,diag,n,A,x,incx);
}
