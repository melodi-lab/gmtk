
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#


if HAVE_HG
# this will keep the Mercurial revision id current if you're building 
# in a working directory. If you're not, use the distributed hgstamp.h
VERSIONFLAG = -DHGID="\"`hg id`\"" -DHGDATE="\"`hg parents --template '{date|date}'`\""
else
VERSIONFLAG = 
endif
.PHONY: testHashDrivers

DEBUGFLAGS =


if CYGWIN_FLAGS
  ANSI = 
else
  ANSI = -ansi -pedantic
endif

# If configure determines we're using GCC C/C++ compilers, set flags appropriately
if USE_GCC_FLAGS
# Only used for GCC compilers
  FEATURE_FLAGS = -D _XOPEN_SOURCE=501
  OPTFLAGS = -O3
  GCC_FLAGS = -Wall $(ANSI)
else
# unknown compiler
  OPTFLAGS =
  GCC_FLAGS =
endif

DBFLAGS =
if NO_ASSERTIONS
# -D is fairly standard accross all C/C++ compilers
DBFLAGS += -DNDEBUG
endif
if NO_MODEL_CHECKS
DBFLAGS += -DNMODELCHECK
endif


# used by clients of miscSupport
noinst_HEADERS = machine-dependent.h

# Build libmiscSupport.a, but don't install it in $(libdir)
noinst_LIBRARIES = libmiscSupport.a

# The source files libmiscSupport.a depends on
libmiscSupport_a_SOURCES = \
hgstamp.h \
counted_ptr.h \
fixed_filter.h \
hash_mtable.h \
matrix.h \
shash_map_iter.h \
strhash_map.h \
error.h error.cc \
popen_err.h popen_err.c \
bp_range.h bp_range.cc \
range.h range.cc \
arguments.h arguments.cc \
logp.h logp.cc \
rand.h rand.cc \
sArray.h sArray.cc \
mArray.h mArray.cc \
cArray.h cArray.cc \
general.h general.cc \
lineqsolve.h lineqsolve.cc \
debug.h debug.cc \
fileParser.h fileParser.cc \
hash_abstract.h hash_abstract.cc \
phipac.h \
mul_mdmd_md.c \
mul_mdmd_md_l0g.c \
mul_mdmd_md_l0nf.c \
mul_mfmf_mf.c \
mul_mfmf_mf_l0g.c \
mul_mfmf_mf_l0nf.c \
eig.h eig.c \
bin_search.cc \
lms_filter.h lms_filter.cc \
rls_filter.h rls_filter.cc \
vhash_set.h \
prime.h prime.c

# GMTKCPPCMD is set to an appropriate default CPP for GMTK
# to process its ASCII input files
libmiscSupport_a_CPPFLAGS = -DDEFAULT_CPP="\"$(GMTKCPPCMD)\"" $(VERSIONFLAG) $(DBFLAGS) $(FEATURE_FLAGS)

# Compilation flags for miscSupport library objects (user's CXXFLAGS will be appended)
libmiscSupport_a_CXXFLAGS = $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS)
libmiscSupport_a_CFLAGS =   $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS)


# Complilation flags for test programs (user's CXXFLAGS will be appended)
AM_CPPFLAGS = $(VERSIONFLAG) $(DBFLAGS) $(FEATURE_FLAGS)
AM_CXXFLAGS = $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS) -DMAIN
AM_CCFLAGS =  $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS) -DMAIN

# Linker flags for test programs (user's LDFLAGS will be appended)
AM_LDFLAGS = 

# Link the test programs with libmiscSupport
LDADD = libmiscSupport.a

# These are test programs, so they only need to be built for make check.
# They are not installed anywhere
check_PROGRAMS = \
testPrime \
testError \
testBPRange \
testPermute \
testPopenErr \
testArguments testLogp testRand testsArray testcArray testmArray \
fileParserTest testHashMapList testHashTree testVHashMap testSHashSet \
testSHashMap testVSHashMap testRLS testLMS testLZERO testQSort
# testVHashSet 

testPermute_SOURCES = testPermute.cc

testQSort_SOURCES = testQSort.cc

# testDebug      doesn't compile?
testError_SOURCES = error.cc error.h
#testDebug_SOURCES = debug.cc debug.h
testBPRange_SOURCES = bp_range.cc bp_range.h

# The sources that comprise each test program
testArguments_SOURCES = arguments.cc arguments.h
testLogp_SOURCES = logp.cc logp.h
testRand_SOURCES = rand.cc rand.h
testsArray_SOURCES = sArray.cc sArray.h
testcArray_SOURCES = cArray.cc cArray.h
testmArray_SOURCES = mArray.cc mArray.h
fileParserTest_SOURCES = fileParser.cc fileParser.h
testLZERO_SOURCES = testLZERO.cc

testHashMapList_SOURCES = hash_map_list.cc hash_map_list.h
testHashTree_SOURCES = hash_tree.h hash_tree.cc
#testVHashSet_SOURCES = vhash_set.h vhash_set.cc
testVHashMap_SOURCES = vhash_map.h vhash_map.cc
testSHashSet_SOURCES = shash_set.h shash_set.cc
testSHashMap_SOURCES = shash_map.h shash_map.cc
testVSHashMap_SOURCES = vshash_map.h vshash_map.cc

testRLS_SOURCES = rls_filter.cc rls_filter.h  adaptive_filter.h
testLMS_SOURCES = lms_filter.cc lms_filter.h  adaptive_filter.h

testHashDrivers: \
testHashMapList testHashTree testVHashSet testVHashMap testSHashSet \
testSHashMap testVSHashMap

