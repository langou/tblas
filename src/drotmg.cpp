#include "blas.h"
#include "rotmg.h"

using tblas::rotmg;

void drotmg_(double &d1, double &d2, double &x1, double &y1, double *param)
{
    int flag=rotmg(d1,d2,x1,y1,param+1);
    param[0]=static_cast<double>(flag);
}
