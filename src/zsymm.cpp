#include "blas.h"
#include "symm.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::symm;

void zsymm_(const char &Side, const char &Uplo, const int &m, const int &n, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC)
{
    const complex<double> zero(0.0);
    const complex<double> one(1.0);
    int info=0;
    char side=toupper(Side);
    char uplo=toupper(Uplo);
    if((side!='L')&&(side!='R'))
        info=1;
    else if((uplo!='U')&&(uplo!='L'))
        info=2;
    else if(m<0)
        info=3;
    else if(n<0)
        info=4;
    else if(ldA<max(1,(side=='L')?m:n))
        info=7;
    else if(ldB<max(1,m))
        info=9;
    else if(ldC<max(1,m))
        info=12;
    if(info>0)
        xerbla_("ZSYMM ",info);
    else
        symm(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
}
