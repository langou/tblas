#include "blas.h"
#include "gerc.h"

using std::conj;
using tblas::gerc;

void cgerc_(const int &m, const int &n, const complex<float> &alpha, complex<float> *x, const int &incx, complex<float> *y, const int &incy, complex<float> *A, const int &ldA)
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
        xerbla_("CGERC ",info);
    else
        gerc(m,n,alpha,x,incx,y,incy,A,ldA);
}
