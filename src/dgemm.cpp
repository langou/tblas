#include "blas.h"
#include "gemm.h"
#include <cctype>
#include <utility>

using std::toupper;
using std::max;
using tblas::gemm;

void dgemm_(const char &TransA, const char &TransB, const int &m, const int &n, const int &k, const double &alpha, double *A, const int &ldA, double *B, const int &ldB, const double &beta, double *C, const int &ldC)
{
    int info=0;
    char transA=toupper(TransA);
    char transB=toupper(TransB);
    if(transA=='C')
        transA='T';
    if(transB=='C')
        transB='T';
    if((transA!='N')&&(transA!='T'))
        info=1;
    else if((transB!='N')&&(transB!='T'))
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
        xerbla_("DGEMM ",info);
    else
        gemm(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
}
