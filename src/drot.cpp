#include "blas.h"
#include "rot.h"

using tblas::rot;

void drot_(const int &n, double *x, const int &incx, double *y, const int &incy, const double &c, const double &s)
{
    if(n>0)
        rot(n,x,incx,y,incy,c,s);
}
