#include "blas.h"
#include "rotm.h"

using tblas::rotm;

void srotm_(const int &n, float *x, const int &incx, float *y, const int &incy, float *param)
{
    if(n>0)
        rotm(n,x,incx,y,incy,param+1,static_cast<int>(param[0]));
}
