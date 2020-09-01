#include "blas.h"
#include "rotm.h"

using tblas::rotm;

void drotm_(const int &n, double *x, const int &incx, double *y, const int &incy, double *param)
{
    if(n>0)
        rotm(n,x,incx,y,incy,param+1,static_cast<int>(param[0]));
}
