#include "blas.h"
#include "rotg.h"

using tblas::rotg;

void zrotg_(complex<double> &a,const complex<double> &b, double &c, complex<double> &s)
{
    rotg(a,b,c,s);
}
