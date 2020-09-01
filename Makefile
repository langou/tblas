####
#### Top level Makefile for TBLAS library
####

CXX=c++
CXXFLAGS=-O3 -w -std=c++11
PREFIX=/usr/local
LIB=libcoblas.a
INC=coblas.h
INSTALL=install

default: blaslib

all: blaslib

blaslib:
	@cd src;make CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)"

install: blaslib
	$(INSTALL) -d $(PREFIX)/lib $(PREFIX)/include
	$(INSTALL) lib/$(LIB) $(PREFIX)/lib
	$(INSTALL) include/$(INC) $(PREFIX)/include

clean:
	@cd src;make clean

