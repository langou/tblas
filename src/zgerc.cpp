#include "blas.h"
#include "gerc.h"

using tblas::gerc;

void zgerc_(const int &m, const int &n, const complex<double> &alpha, complex<double> *x, const int &incx, complex<double> *y, const int &incy, complex<double> *A, const int &ldA)
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
        xerbla_("ZGERC ",info);
    else
        gerc(m,n,alpha,x,incx,y,incy,A,ldA);
}
