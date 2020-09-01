#include "blas.h"
#include "herk.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::herk;

void zherk_(const char &Uplo, const char &Trans, const int &n, const int &k, const double &alpha, complex<double> *A, const int &ldA, const double &beta, complex<double> *C, const int &ldC)
{
    int info=0;
    char trans=toupper(Trans);
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if((trans!='N')&&(trans!='C'))
        info=2;
    else if(n<0)
        info=3;
    else if(k<0)
        info=4;
    else if(ldA<max(1,(trans=='N')?n:k))
        info=7;
    else if(ldC<max(1,n))
        info=10;
    if(info>0)
        xerbla_("ZHERK ",info);
    else
        herk(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
}
