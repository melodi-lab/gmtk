
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2013 Jeff Bilmes
# Licensed under the Open Software License version 3.0
#


DEBUGFLAGS = -ggdb3
#OPTFLAGS = -O3 -march=native -mtune=native -mfpmath=sse -msse3 -flto
OPTFLAGS = -O3 -march=native -mtune=native -mfpmath=sse -msse3
ANSI = -std=c++0x
#PEDANTIC = -pedantic
#WARNING = -Wall

AM_CPPFLAGS = $(BLAS_CPPFLAGS) $(BOOST_CPPFLAGS) -D__LINUX

AM_CXXFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(PEDANTIC) $(WARNING) $(ANSI)

#AM_LDFLAGS = -static -fwhole-program -L. 
AM_LDFLAGS = -L. $(BOOST_LDFLAGS) $(BLAS_LDFLAGS)

LDADD = libDMLP.a

LIBS += -lboost_thread -lboost_filesystem -lboost_system $(BLAS_LIBS)
# -lpthread

noinst_LIBRARIES = libDMLP.a
noinst_bin_PROGRAMS = dmlp

dmlp_SOURCES = dmlp.cc

libDMLP_a_SOURCES =           \
  Globals.h                   \
  DBN.h Layer.h TrainSGD.h    \
  Matrix.h Matrix.cpp         \
  MatrixFunc.h MatrixFunc.cc  \
  Random.h Random.cc
