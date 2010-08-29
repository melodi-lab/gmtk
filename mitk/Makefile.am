
AM_CPPFLAGS = \
-I$(srcdir)/../miscSupport \
-I$(srcdir)/../featureFileIO \
-I$(srcdir)/../IEEEFloatingpoint

LDADD = \
$(builddir)/../featureFileIO/libgmtkio.a \
$(builddir)/../miscSupport/libmiscSupport.a \
$(builddir)/../IEEEFloatingpoint/libIEEEsupport.a

EXCCFLAGS = 
EXCXXFLAGS = 
EXLDFLAGS  = 
DEBUGFLAGS = -g
OPTFLAGS = -O3

ANSI = -ansi
PEDANTIC = -pedantic -Wno-long-long
WALL = -Wall -Wno-deprecated

AM_CFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)
AM_CXXFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

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
