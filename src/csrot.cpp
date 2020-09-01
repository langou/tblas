#include "blas.h"
#include "rot.h"

using tblas::rot;

void csrot_(const int &n, complex<float> *x, const int &incx, complex<float> *y, const int &incy, const float &c, const float &s)
{
    if(n>0)
        rot(n,x,incx,y,incy,c,s);
}
