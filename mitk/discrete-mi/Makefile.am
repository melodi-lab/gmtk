
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
-I$(srcdir)/../../miscSupport \
-I$(srcdir)/../../featureFileIO \
-I$(builddir)/../../IEEEFloatingpoint \
$(DBFLAGS)

LDADD = \
$(builddir)/../../miscSupport/libmiscSupport.a \
$(builddir)/../../featureFileIO/libgmtkio.a

AM_CXXFLAGS = $(DEBUGFLAGS) $(OPTFLAGS) $(GCC_FLAGS)

bin_PROGRAMS = discrete-mi
discrete_mi_SOURCES = \
mixNormalCollection.h mixNormalCollection.cc \
discrete_mi.cc \
mixNormal.h mixNormal.cc \
readRange.h readRange.cc \
rangeSetCollection.cc

