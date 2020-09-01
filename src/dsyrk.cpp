#include "blas.h"
#include "syrk.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::syrk;

void dsyrk_(const char &Uplo, const char &Trans, const int &n, const int &k, const double &alpha, double *A, const int &ldA, const double &beta, double *C, const int &ldC)
{
    int info=0;
    char trans=toupper(Trans);
    char uplo=toupper(Uplo);
    if(trans=='C')
        trans='T';
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
    else if(ldC<max(1,n))
        info=10;
    if(info>0)
        xerbla_("DSYRK ",info);
    else
        syrk(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
}
