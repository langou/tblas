#include "blas.h"
#include "rot.h"

using tblas::rot;

void zdrot_(const int &n, complex<double> *x, const int &incx, complex<double> *y, const int &incy, const double &c, const double &s)
{
    if(n>0)
        rot(n,x,incx,y,incy,c,s);
}
