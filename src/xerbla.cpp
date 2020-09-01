#include "blas.h"
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;
using std::exit;

void xerbla_(const char *name, const int &info)
{
    cout << " ** On entry to " << name;
    cout << " parameter number " << info;
    cout << " had an illegal value." << endl;
    exit(EXIT_FAILURE);
}
