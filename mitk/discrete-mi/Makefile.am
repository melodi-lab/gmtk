
AM_CPPFLAGS = \
-I$(srcdir)/../../miscSupport \
-I$(srcdir)/../../featureFileIO \
-I$(srcdir)/../../IEEEFloatingpoint

LDADD = \
$(builddir)/../../miscSupport/libmiscSupport.a \
$(builddir)/../../featureFileIO/libgmtkio.a


EXCCFLAGS = 
EXCXXFLAGS = 
EXLDFLAGS  = 
DEBUGFLAGS = -g
OPTFLAGS = -O3

ANSI = 
PEDANTIC = 
WALL = 

AM_CXXFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

EXTRA_PROGRAMS = discrete-mi
discrete_mi_SOURCES = \
mixNormalCollection.h mixNormalCollection.cc \
discrete_mi.h discrete_mi.cc \
mixNormal.h mixNormal.cc \
readRange.h readRange.cc \
rangeSetCollection.h rangeSetCollection.cc








