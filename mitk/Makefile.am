
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#


DEBUGFLAGS = 


# If configure determines we're using GCC C/C++ compilers, set flags appropriately
if USE_GCC_FLAGS
# Only used for GCC compilers
  OPTFLAGS = -O3
  GCC_FLAGS = -Wall -ansi -pedantic
else
# unknown compiler
  OPTFLAGS =
  GCC_FLAGS =
endif

DBFLAGS =
if NO_ASSERTIONS
DBFLAGS += -DNDEBUG
endif
if NO_MODEL_CHECKS
DBFLAGS += -DNMODELCHECK
endif


AM_CPPFLAGS = \
-I$(srcdir)/../miscSupport \
-I$(srcdir)/../featureFileIO \
-I$(builddir)/../IEEEFloatingpoint \
$(DBFLAGS)

LDADD = \
$(builddir)/../featureFileIO/libgmtkio.a \
$(builddir)/../miscSupport/libmiscSupport.a \
$(builddir)/../IEEEFloatingpoint/libIEEEsupport.a

AM_CFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(GCC_FLAGS)
AM_CXXFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(GCC_FLAGS)

SUBDIRS = discrete-mi

EXTRA_DIST = global-parameters.h mixNormal.h
noinst_PROGRAMS = multivariate-mi bivariate-mi
check_PROGRAMS = mvnrnd testReadRange test-matrix-ops

mvnrnd_SOURCES = mvnrnd.cc matrix-ops.h matrix-ops.cc

testReadRange_SOURCES = testReadRange.cc \
readRange.h readRange.cc \
rangeSetCollection.cc

multivariate_mi_SOURCES = multivariate_mi.cc \
mixNormalCollection.h mixNormalCollection.cc \
mixNormal.cc \
kmeans.cc \
readRange.h readRange.cc \
rangeSetCollection.cc \
partition.cc \
computeMI.cc \
computeMI-data.cc \
data-points.h data-points.cc \
matrix-ops.h matrix-ops.cc \
tests.h tests.cc \
io.cc \
generate-synthetic-data.cc \
compute-entropy.cc

bivariate_mi_SOURCES = bivariate_mi.cc \
MixBiNormal.h MixBiNormal.cc \
MixBiNormal_chelp.c

test_matrix_ops_SOURCES = test-matrix-ops.cc \
matrix-ops.h matrix-ops.cc

