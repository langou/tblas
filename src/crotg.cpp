#include "blas.h"
#include "rotg.h"

using tblas::rotg;

void crotg_(complex<float> &a, const complex<float> &b, float &c, complex<float> &s)
{
    rotg(a,b,c,s);
}
