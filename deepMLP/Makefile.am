
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2013 Jeff Bilmes
# Licensed under the Open Software License version 3.0
#


DEBUGFLAGS = -ggdb3
#OPTFLAGS = -O3 -march=native -mtune=native -mfpmath=sse -msse3 -flto
OPTFLAGS = -O3 -march=native -mtune=native -mfpmath=sse -msse3
ANSI = 
#PEDANTIC = -pedantic
#WARNING = -Wall

AM_CPPFLAGS = \
-I$(srcdir)/../miscSupport   \
-I$(srcdir)/../featureFileIO

AM_CXXFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(PEDANTIC) $(WARNING) $(ANSI) $(BLAS_CFLAGS)

#AM_LDFLAGS = -static -fwhole-program -L. 
AM_LDFLAGS = -L.

LDADD =                                     \
libDMLP.a                                   \
$(builddir)/../featureFileIO/libgmtkio.a    \
$(builddir)/../miscSupport/libmiscSupport.a 

if BUILD_PHIPAC
LDADD += libPHiPAC.a
endif

LIBS += $(BLAS_LIBS)

noinst_LIBRARIES = $(dmlp_library) $(optional_libraries)
EXTRA_LIBRARIES  = libDMLP.a libPHiPAC.a

EXTRA_PROGRAMS  = dmlp mnist label dmlpTest orderTest phipacTest

check_PROGRAMS = $(check_programs)
#dmlpTest orderTest phipacTest flat2vit
dmlp_SOURCES = dmlp.cc MNIST.h
mnist_SOURCES = mnist.cc MNIST.h
label_SOURCES = label.cc
dmlpTest_SOURCES = dmlpTest.cc
orderTest_SOURCES = orderTest.cc
phipacTest_SOURCES = phipacTest.cc


libDMLP_a_SOURCES =           \
  dummy.cc Globals.h          \
  miniblas.h miniblas.cc      \
  DBN.h DBN.cc Layer.h        \
  Matrix.h Matrix.cc          \
  MatrixFunc.h MatrixFunc.cc  \
  MMapMatrix.h MMapMatrix.cc  \
  BatchSource.h               \
  AsynchronousBatchSource.h   \
  AsynchronousBatchSource.cc


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
