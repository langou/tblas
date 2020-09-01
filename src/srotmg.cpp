#include "blas.h"
#include "rotmg.h"

using tblas::rotmg;

void srotmg_(float &d1, float &d2, float &x1, float &y1, float *param)
{
    int flag=rotmg(d1,d2,x1,y1,param+1);
    param[0]=static_cast<float>(flag);
}
