#include "blas.h"
#include "gemm.h"
#include <cctype>
#include <utility>

using std::max;
using std::toupper;
using tblas::gemm;

void sgemm_(const char &TransA, const char &TransB, const int &m, const int &n, const int &k, const float &alpha, float *A, const int &ldA, float *B, const int &ldB, const float &beta, float *C, const int &ldC)
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
        xerbla_("SGEMM ",info);
    else
        gemm(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
}
