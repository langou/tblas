#include "blas.h"
#include "spr2.h"
#include <cctype>

using std::toupper;
using tblas::spr2;

void dspr2_(const char &Uplo, const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy, double *A)
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
    if(info>0)
        xerbla_("DSPR2 ",info);
    else
        spr2(uplo,n,alpha,x,incx,y,incy,A);
}
