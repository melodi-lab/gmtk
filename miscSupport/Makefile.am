
if HAVE_HG
# this will keep the Mercurial revision id current if you're building 
# in a working directory. If you're not, use the distributed hgstamp.h
VERSIONFLAG = -DHGID="\"`hg id`\"" -DHGDATE="\"`hg parents --template '{date|date}'`\""
else
VERSIONFLAG = 
endif
.PHONY: testHashDrivers

EXCCFLAGS = 
EXCXXFLAGS = 
EXLDFLAGS  = 
DEBUGFLAGS = -g
if NO_ASSERTIONS
DEBUGFLAGS += -DNDEBUG
endif
if NO_MODEL_CHECKS
DEBUGFLAGS += -DNMODELCHECK
endif
OPTFLAGS = -O3
HASH = 

if CYGWIN_FLAGS
ANSI = 
PEDANTIC =
else
ANSI = -ansi
PEDANTIC = -pedantic
endif

WALL = -Wall

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
mul_mdmd_md.c \
mul_mdmd_md_l0g.c \
mul_mdmd_md_l0nf.c \
eig.h eig.c \
bin_search.cc \
lms_filter.h lms_filter.cc \
rls_filter.h rls_filter.cc \
vhash_set.h

libmiscSupport_a_CPPFLAGS = -DDEFAULT_CPP="\"$(GMTKCPPCMD)\""

# Compilation flags for miscSupport library objects (user's CXXFLAGS will
# be appended)
libmiscSupport_a_CXXFLAGS = $(VERSIONFLAG) \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC) $(HASH)
libmiscSupport_a_CFLAGS = $(VERSIONFLAG) \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCCFLAGS) $(WALL) $(ANSI) $(PEDANTIC) $(HASH)


# Complilation flags for test programs (user's CXXFLAGS will be appended)
AM_CXXFLAGS = $(VERSIONFLAG) $(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC) $(HASH) -DMAIN

# Linker flags for test programs (user's LDFLAGS will be appended)
AM_LDFLAGS = $(EXLDFLAGS)

# Link the test programs with libmiscSupport
LDADD = libmiscSupport.a

# These are test programs, so they only need to be built for make check.
# They are not installed anywhere
check_PROGRAMS = \
testError \
testBPRange \
testArguments testLogp testRand testsArray testcArray testmArray \
fileParserTest testHashMapList testHashTree testVHashMap testSHashSet \
testSHashMap testVSHashMap testRLS testLMS testLZERO
# testVHashSet 

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

