#include "blas.h"
#include "syr2k.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::syr2k;

void zsyr2k_(const char &Uplo, const char &Trans, const int &n, const int &k, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC)
{
    int info=0;
    char trans=toupper(Trans);
    char uplo=toupper(Uplo);
    if((uplo!='U')&&(uplo!='L'))
        info=1;
    else if((trans!='N')&&(trans!='T'))
        info=2;
    else if(n<0)
        info=3;
    else if(k<0)
        info=4;
    else if(ldA<max(1,(trans=='N')?n:k))
        info=7;
    else if(ldB<max(1,(trans=='N')?n:k))
        info=9;
    else if(ldC<max(1,n))
        info=12;
    if(info>0)
        xerbla_("ZSYR2K ",info);
    else
        syr2k(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
}
