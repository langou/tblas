####
#### Makefile for legacy interface to TBLAS
####

CXX=c++
CXXFLAGS=-O3 -w -std=c++11
TARGET=libtblas.a
LIBTOOL=ar cr
RANLIB=ranlib
INSTALL=install
INCDIR=../include
LIBDIR=../lib
LIB=

OBJ=caxpy.o ccopy.o cdotc.o cdotu.o cgbmv.o cgemm.o cgemv.o cgerc.o cgeru.o chbmv.o \
chemm.o chemv.o cher.o cher2.o cher2k.o cherk.o chpmv.o chpr.o chpr2.o crotg.o \
cscal.o csrot.o csscal.o cswap.o csymm.o csyr2k.o csyrk.o ctbmv.o ctbsv.o ctpmv.o \
ctpsv.o ctrmm.o ctrmv.o ctrsm.o ctrsv.o dasum.o daxpy.o dcopy.o ddot.o \
dgbmv.o dgemm.o dgemv.o dger.o dnrm2.o drot.o drotg.o drotm.o drotmg.o dsbmv.o \
dscal.o dsdot.o dspmv.o dspr.o dspr2.o dswap.o dsymm.o dsymv.o dsyr.o dsyr2.o \
dsyr2k.o dsyrk.o dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o dtrmv.o dtrsm.o \
dtrsv.o dzasum.o dznrm2.o icamax.o idamax.o isamax.o izamax.o sasum.o saxpy.o \
scasum.o scnrm2.o scopy.o sdot.o sdsdot.o sgbmv.o sgemm.o sgemv.o sger.o snrm2.o \
srot.o srotg.o srotm.o srotmg.o ssbmv.o sscal.o sspmv.o sspr.o sspr2.o sswap.o \
ssymm.o ssymv.o ssyr.o ssyr2.o ssyr2k.o ssyrk.o stbmv.o stbsv.o stpmv.o stpsv.o \
strmm.o strmv.o strsm.o strsv.o zaxpy.o zcopy.o zdotc.o zdotu.o \
zdrot.o zdscal.o zgbmv.o zgemm.o zgemv.o zgerc.o zgeru.o zhbmv.o zhemm.o zhemv.o \
zher.o zher2.o zher2k.o zherk.o zhpmv.o zhpr.o zhpr2.o zrotg.o zscal.o zswap.o \
zsymm.o zsyr2k.o zsyrk.o ztbmv.o ztbsv.o ztpmv.o ztpsv.o ztrmm.o ztrmv.o \
ztrsm.o ztrsv.o xerbla.o

default: $(TARGET)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(INSTALL) -d $(LIBDIR)
	$(LIBTOOL) $(LIBDIR)/$(TARGET) $(OBJ) $(LIB)
	$(RANLIB) $(LIBDIR)/$(TARGET)

saxpy.o caxpy.o daxpy.o zaxpy.o: $(INCDIR)/axpy.h
sasum.o scasum.o dasum.o dzasum.o: $(INCDIR)/asum.h
scopy.o ccopy.o dcopy.o zcopy.o: $(INCDIR)/copy.h
sdot.o cdotu.o ddot.o zdotu.o sdsdot.o dsdot.o: $(INCDIR)/dot.h
cdotc.o zdotc.o: $(INCDIR)/dotc.h
sgbmv.o cgbmv.o dgbmv.o zgbmv.o: $(INCDIR)/gbmv.h
sgemm.o cgemm.o dgemm.o zgemm.o: $(INCDIR)/gemm.h
sgemv.o cgemv.o dgemv.o zgemv.o: $(INCDIR)/gemv.h
sger.o cgeru.o dger.o zgeru.o: $(INCDIR)/ger.h
cgerc.o zgerc.o: $(INCDIR)/gerc.h
chbmv.o zhbmv.o: $(INCDIR)/hbmv.h
chemm.o zhemm.o: $(INCDIR)/hemm.h
chemv.o zhemv.o: $(INCDIR)/hemv.h
cher.o zher.o: $(INCDIR)/her.h
cher2.o zher2.o: $(INCDIR)/her2.h
cher2k.o zher2k.o: $(INCDIR)/her2k.h
cherk.o zherk.o: $(INCDIR)/herk.h
chpmv.o zhpmv.o: $(INCDIR)/hpmv.h
chpr.o zhpr.o: $(INCDIR)/hpr.h
chpr2.o zhpr2.o: $(INCDIR)/hpr2.h
isamax.o icamax.o idamax.o izamax.o: $(INCDIR)/imax.h
snrm2.o dnrm2.o scnrm2.o dznrm2.o: $(INCDIR)/nrm2.h
srot.o drot.o csrot.o zdrot.o: $(INCDIR)/rot.h
srotg.o crotg.o drotg.o zrotg.o: $(INCDIR)/rotg.h
ssbmv.o dsbmv.o: $(INCDIR)/sbmv.h
sscal.o csscal.o cscal.o dscal.o zdscal.o zscal.o: $(INCDIR)/scal.h
sspmv.o dspmv.o: $(INCDIR)/spmv.h
sspr.o dspr.o: $(INCDIR)/spr.h
sspr2.o dspr2.o: $(INCDIR)/spr2.h
sswap.o cswap.o dswap.o zswap.o: $(INCDIR)/swap.h
ssymm.o csymm.o dsymm.o zsymm.o: $(INCDIR)/symm.h
ssymv.o dsymv.o: $(INCDIR)/symv.h
ssyr.o dsyr.o: $(INCDIR)/syr.h
ssyr2.o dsyr2.o: $(INCDIR)/syr2.h
ssyr2k.o csyr2k.o dsyr2k.o zsyr2k.o: $(INCDIR)/syr2k.h
ssyrk.o csyrk.o dsyrk.o zsyrk.o: $(INCDIR)/syrk.h
stbmv.o ctbmv.o dtbmv.o ztbmv.o: $(INCDIR)/tbmv.h
stbsv.o ctbsv.o dtbsv.o ztbsv.o: $(INCDIR)/tbsv.h
stpmv.o ctpmv.o dtpmv.o ztpmv.o: $(INCDIR)/tpmv.h
stpsv.o ctpsv.o dtpsv.o ztpsv.o: $(INCDIR)/tpsv.h
strmm.o ctrmm.o dtrmm.o ztrmm.o: $(INCDIR)/trmm.h
strmv.o ctrmv.o dtrmv.o ztrmv.o: $(INCDIR)/trmv.h
strsm.o ctrsm.o dtrsm.o ztrsm.o: $(INCDIR)/trsm.h
strsv.o ctrsv.o dtrsv.o ztrsv.o: $(INCDIR)/trsv.h

.cpp.o:
	$(CXX) -c $(CXXFLAGS) -I$(INCDIR) $<

clean:
	rm -f $(OBJ) $(LIBDIR)/$(TARGET)
