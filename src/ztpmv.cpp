#include "blas.h"
#include "tpmv.h"
#include <cctype>

using std::toupper;
using tblas::tpmv;

void ztpmv_(const char &Uplo, const char &Trans, const char &Diag, const int &n, complex<double> *A, complex<double> *x, const int &incx)
{
    int info=0;
    char uplo=toupper(Uplo);
    char trans=toupper(Trans);
    char diag=toupper(Diag);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if((trans!='N')&&(trans!='T')&&(trans!='C'))
        info=2;
    else if((diag!='U')&&(diag!='N'))
        info=3;
    else if(n<0)
        info=4;
    else if(incx==0)
        info=7;
    if(info>0)
        xerbla_("ZTPMV ",info);
    else
        tpmv(uplo,trans,diag,n,A,x,incx);
}
