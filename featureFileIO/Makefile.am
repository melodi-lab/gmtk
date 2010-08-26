
EXCCFLAGS = 
EXCXXFLAGS = 
EXLDFLAGS = 
DEBUGFLAGS = -g
OPTFLAGS = -O3

ANSI = -ansi
PEDANTIC = -pedantic
WALL = -Wall


AM_CPPFLAGS = -I$(srcdir)/../miscSupport -I$(srcdir)/../IEEEFloatingpoint

AM_CFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCCFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

AM_CXXFLAGS = \
$(DEBUGFLAGS) $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

AM_LDFLAGS = $(EXLDFLAGS)

LDADD = \
libgmtkio.a \
$(srcdir)/../miscSupport/libmiscSupport.a \
$(srcdir)/../IEEEFloatingpoint/libIEEEsupport.a

noinst_LIBRARIES = libgmtkio.a

libgmtkio_a_SOURCES = \
GMTK_ObservationMatrix.h GMTK_ObservationMatrix.cc \
GMTK_Stream.h GMTK_Stream.cc \
GMTK_DataTransformations.h GMTK_DataTransformations.cc \
GMTK_WordOrganization.h GMTK_WordOrganization.cc \
GMTK_ObsStats.h GMTK_ObsStats.cc \
pfile.h pfile.cc \
vbyteswapping.h vbyteswapping.cc

libgmtkio_a_LIBADD = \
$(srcdir)/../miscSupport/libmiscSupport.a \
$(srcdir)/../IEEEFloatingpoint/libIEEEsupport.a

noinst_PROGRAMS = obs-print obs-info obs-skmeans obs-concat obs-diff
check_PROGRAMS = gmtest

gmtest_SOURCES = gmtest.cc

obs_print_SOURCES = \
GMTK_ObsPrint.h GMTK_ObsPrint.cc \
GMTK_ObsKLT.h GMTK_ObsKLT.cc \
GMTK_ObsNorm.h GMTK_ObsNorm.cc \
GMTK_ObsGaussianNorm.h GMTK_ObsGaussianNorm.cc \
GMTK_ObsAddSil.h GMTK_ObsAddSil.cc

obs_info_SOURCES = GMTK_ObsInfo.h  GMTK_ObsInfo.cc

obs_skmeans_SOURCES = \
GMTK_ObsSkmeans.h GMTK_ObsSkmeans.cc \
GMTK_ObsInitmg.h GMTK_ObsInitmg.cc \
GMTK_Kmeans.h GMTK_Kmeans.cc

obs_concat_SOURCES = GMTK_ObsConcat.h GMTK_ObsConcat.cc

obs_diff_SOURCES = GMTK_ObsDiff.h GMTK_ObsDiff.cc


