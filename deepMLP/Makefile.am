
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

AM_CPPFLAGS = $(BLAS_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(srcdir)/../miscSupport

AM_CXXFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(PEDANTIC) $(WARNING) $(ANSI)

#AM_LDFLAGS = -static -fwhole-program -L. 
AM_LDFLAGS = -L. $(BOOST_LDFLAGS) $(BLAS_LDFLAGS)

if USE_PHIPAC
LDADD = libDMLP.a libPHiPAC.a
else
LDADD = libDMLP.a
endif

LIBS += -lboost_thread -lboost_filesystem -lboost_system $(BLAS_LIBS)
# -lpthread

noinst_LIBRARIES = libDMLP.a libPHiPAC.a

noinst_PROGRAMS = dmlp mnist label

dmlp_SOURCES = dmlp.cc MNIST.h

mnist_SOURCES = mnist.cc MNIST.h

label_SOURCES = label.cc

check_PROGRAMS = dmlpTest orderTest phipacTest

dmlpTest_SOURCES = dmlpTest.cc
orderTest_SOURCES = orderTest.cc
phipacTest_SOURCES = phipacTest.cc


libDMLP_a_SOURCES =           \
  dummy.cc Globals.h          \
  miniblas.h miniblas.cc      \
  DBN.h Layer.h               \
  Matrix.h Matrix.cc          \
  MatrixFunc.h MatrixFunc.cc  \
  Random.h Random.cc

libPHiPAC_a_SOURCES =     \
mm_double_NN_1.c          \
mm_double_NN_1_general.c  \
mm_double_NN_1_l1.c       \
mm_double_NN_1_nofringe.c \
mm_double_NN_c.c          \
mm_double_NN_c_general.c  \
mm_double_NN_c_l1.c       \
mm_double_NN_c_nofringe.c \
mm_double_NT_1.c          \
mm_double_NT_1_general.c  \
mm_double_NT_1_l1.c       \
mm_double_NT_1_nofringe.c \
mm_double_NT_c.c          \
mm_double_NT_c_general.c  \
mm_double_NT_c_l1.c       \
mm_double_NT_c_nofringe.c \
mm_double_TN_1.c          \
mm_double_TN_1_general.c  \
mm_double_TN_1_l1.c       \
mm_double_TN_1_nofringe.c \
mm_double_TN_c.c          \
mm_double_TN_c_general.c  \
mm_double_TN_c_l1.c       \
mm_double_TN_c_nofringe.c \
mm_double_TT_1.c          \
mm_double_TT_1_general.c  \
mm_double_TT_1_l1.c       \
mm_double_TT_1_nofringe.c \
mm_double_TT_c.c          \
mm_double_TT_c_general.c  \
mm_double_TT_c_l1.c       \
mm_double_TT_c_nofringe.c \
phipac_dgemm.c 
