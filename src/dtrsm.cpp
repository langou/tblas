#include "blas.h"
#include "trsm.h"
#include <cctype>
#include <utility>

using std::toupper;
using std::max;
using tblas::trsm;

void dtrsm_(const char &Side, const char &Uplo, const char &Trans, const char &Diag, const int &m, const int &n, const double &alpha, double *A, const int &ldA, double *B, const int &ldB)
{
    int info=0;
    char side=toupper(Side);
    char uplo=toupper(Uplo);
    char trans=toupper(Trans);
    char diag=toupper(Diag);
    if(trans=='C')
        trans='T';
    if((side!='L')&&(side!='R'))
        info=1;
    else if((uplo!='U')&&(uplo!='L'))
        info=2;
    else if((trans!='N')&&(trans!='T')&&(trans!='C'))
        info=3;
    else if((diag!='U')&&(diag!='N'))
        info=4;
    else if(m<0)
        info=5;
    else if(n<0)
        info=6;
    else if(ldA<max(1,(side=='L')?m:n))
        info=9;
    else if(ldB<max(1,m))
        info=11;
    if(info!=0)
        xerbla_("DTRSM ",info);
    else
        trsm(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
}
