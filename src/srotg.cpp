#include "blas.h"
#include "rotg.h"

using tblas::rotg;

void srotg_(float &a, float &b, float &c, float &s)
{
    rotg(a,b,c,s);
}

