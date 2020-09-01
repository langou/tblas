#include "blas.h"
#include "her.h"
#include <cctype>

using std::toupper;
using tblas::her;

void cher_(const char &Uplo, const int &n, const float &alpha, complex<float> *x, const int &incx, complex<float> *A, const int &ldA)
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
        xerbla_("CHER  ",info);
    else
        her(uplo,n,alpha,x,incx,A,ldA);
}
