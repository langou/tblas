#include "blas.h"
#include "her2k.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::her2k;

void cher2k_(const char &Uplo, const char &Trans, const int &n, const int &k, const complex<float> &alpha, complex<float> *A, const int &ldA, complex<float> *B, const int &ldB, const float &beta, complex<float> *C, const int &ldC)
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
    else if(ldB<max(1,(trans=='N')?n:k))
        info=9;
    else if(ldC<max(1,n))
        info=12;
    if(info>0)
        xerbla_("CHER2K ",info);
    else
        her2k(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
}
