#include "blas.h"
#include "rot.h"

using tblas::rot;

void srot_(const int &n, float *x, const int &incx, float *y, const int &incy, const float &c, const float &s)
{
    if(n>0)
        rot(n,x,incx,y,incy,c,s);
}
