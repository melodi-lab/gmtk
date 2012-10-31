
AM_CPPFLAGS = \
-I$(srcdir)/../../miscSupport \
-I$(srcdir)/../../featureFileIO \
-I$(builddir)/../../IEEEFloatingpoint

LDADD = \
$(builddir)/../../miscSupport/libmiscSupport.a \
$(builddir)/../../featureFileIO/libgmtkio.a


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

ANSI = -ansi 
PEDANTIC = -pedantic
WALL = -Wall

AM_CXXFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

EXTRA_PROGRAMS = discrete-mi
discrete_mi_SOURCES = \
mixNormalCollection.h mixNormalCollection.cc \
discrete_mi.cc \
mixNormal.h mixNormal.cc \
readRange.h readRange.cc \
rangeSetCollection.cc

