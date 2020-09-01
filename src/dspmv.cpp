#include "blas.h"
#include "spmv.h"
#include <cctype>

using std::toupper;
using tblas::spmv;

void dspmv_(const char &Uplo, const int &n, const double &alpha, double *A, double *x, const int &incx, const double &beta, double *y, const int &incy)
{
    int info=0;
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=6;
    else if(incy==0)
        info=9;
    if(info>0)
        xerbla_("DSPMV ",info);
    else
        spmv(uplo,n,alpha,A,x,incx,beta,y,incy);
}
