#include "blas.h"
#include "gemm.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::gemm;

void zgemm_(const char &TransA, const char &TransB, const int &m, const int &n, const int &k, const complex<double> &alpha, complex<double> *A, const int &ldA, complex<double> *B, const int &ldB, const complex<double> &beta, complex<double> *C, const int &ldC)
{
    int info=0;
    char transA=toupper(TransA);
    char transB=toupper(TransB);
    if((transA!='N')&&(transA!='T')&&(transA!='C'))
        info=1;
    else if((transB!='N')&&(transB!='T')&&(transB!='C'))
        info=2;
    else if(m<0)
        info=3;
    else if(n<0)
        info=4;
    else if(k<0)
        info=5;
    else if(ldA<max(1,(transA!='N')?k:m))
        info=8;
    else if(ldB<max(1,(transB!='N')?n:k))
        info=10;
    else if(ldC<max(1,m))
        info=13;
    if(info>0)
        xerbla_("ZGEMM ",info);
    else
        gemm(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
}
