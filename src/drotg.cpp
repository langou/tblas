#include "blas.h"
#include "rotg.h"

using tblas::rotg;

void drotg_(double &a, double &b, double &c, double &s)
{
    rotg(a,b,c,s);
}
