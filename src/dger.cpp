#include "blas.h"
#include "ger.h"

using tblas::ger;

void dger_(const int &m, const int &n, const double &alpha, double *x, const int &incx, double *y, const int &incy, double *A, const int &ldA)
{
    int info=0;
    if(m<0)
        info=1;
    else if(n<0)
        info=2;
    else if(incx==0)
        info=5;
    else if(incy==0)
        info=7;
    else if(ldA<m)
        info=9;
    if(info!=0)
        xerbla_("DGER  ",info);
    else
        ger(m,n,alpha,x,incx,y,incy,A,ldA);
}
