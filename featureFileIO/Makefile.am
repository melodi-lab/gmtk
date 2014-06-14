
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

DEBUGFLAGS = 


if CYGWIN_FLAGS
  ANSI = 
else
  ANSI = -ansi -pedantic
endif

# If configure determines we're using GCC C/C++ compilers, set flags appropriately
if USE_GCC_FLAGS
# Only used for GCC compilers
  OPTFLAGS = -O3
  GCC_FLAGS = -Wall $(ANSI)
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


AM_CPPFLAGS = -I$(srcdir)/../miscSupport -I$(builddir)/../IEEEFloatingpoint $(VERSIONFLAG) $(DBFLAGS)

AM_CFLAGS = $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS)

AM_CXXFLAGS = $(OPTFLAGS) $(DEBUGFLAGS) $(GCC_FLAGS)

AM_LDFLAGS =

LDADD = \
libgmtkio.a \
$(builddir)/../miscSupport/libmiscSupport.a \
$(builddir)/../IEEEFloatingpoint/libIEEEsupport.a

noinst_LIBRARIES = libgmtkio.a

libgmtkio_a_SOURCES = \
GMTK_ObservationArguments.h ObsArguments.h \
GMTK_ObservationMatrix.h GMTK_ObservationMatrix.cc \
GMTK_Stream.h GMTK_Stream.cc \
GMTK_DataTransformations.cc \
GMTK_WordOrganization.h GMTK_WordOrganization.cc \
GMTK_ObsStats.h GMTK_ObsStats.cc \
pfile.h pfile.cc \
vbyteswapping.h vbyteswapping.cc \
file_utils.h file_utils.cc \
GMTK_ObservationSource.h \
GMTK_StreamSource.h GMTK_StreamSource.cc \
GMTK_ObservationStream.h \
GMTK_ASCIIStream.h GMTK_ASCIIStream.cc \
GMTK_BinStream.h GMTK_BinStream.cc \
GMTK_FilterStream.h \
GMTK_MergeStream.h \
GMTK_FileStream.h GMTK_FileSrcStream.h \
GMTK_FileSource.h GMTK_FileSource.cc \
GMTK_FileSourceNoCache.h GMTK_FileSourceNoCache.cc \
GMTK_CreateFileSource.h GMTK_CreateFileSource.cc \
GMTK_ObservationFile.h GMTK_ObservationFile.cc \
GMTK_PFileFile.h GMTK_PFileFile.cc \
GMTK_HTKFile.h GMTK_HTKFile.cc \
GMTK_HDF5File.h GMTK_HDF5File.cc \
GMTK_ASCIIFile.h GMTK_ASCIIFile.cc \
GMTK_FlatASCIIFile.h GMTK_FlatASCIIFile.cc \
GMTK_BinaryFile.h GMTK_BinaryFile.cc \
GMTK_Filter.h GMTK_Filter.cc GMTK_SubmatrixDescriptor.h \
GMTK_MergeFile.h GMTK_MergeFile.cc \
GMTK_FilterFile.h \
GMTK_FIRFilter.h GMTK_FIRFilter.cc \
GMTK_ARMAFilter.h GMTK_ARMAFilter.cc \
GMTK_UpsampleFilter.h GMTK_UpsampleFilter.cc \
GMTK_UpsampleSmoothFilter.h GMTK_UpsampleSmoothFilter.cc \
GMTK_AffineFilter.h GMTK_AffineFilter.cc \
GMTK_TrainingSchedule.h GMTK_TrainingSchedule.cc \
GMTK_RandomSampleSchedule.h \
GMTK_LinearSchedule.h \
GMTK_PermutationSchedule.h \
GMTK_ShuffleSchedule.h

check_PROGRAMS = \
testPFileFile \
testASCIIFile \
testASCIIStream \
testBinStream \
testFileStream \
testObsFile \
testFileSrc \
testFileWriter \
testPFileWriter \
testFilter \
testStreamSource \
testMerge

#testFileSource 


testPFileFile_SOURCES = testPFileFile.cc 
testASCIIFile_SOURCES = testASCIIFile.cc 
#testFileSource_SOURCES = testFileSource.cc
testASCIIStream_SOURCES = testASCIIStream.cc
testBinStream_SOURCES = testBinStream.cc
testFileStream_SOURCES = testFileStream.cc
testFileSrc_SOURCES = testFileSrc.cc
testFileWriter_SOURCES = testFileWriter.cc
testPFileWriter_SOURCES = testPFileWriter.cc
testFilter_SOURCES = testFilter.cc
testObsFile_SOURCES = testObsFile.cc
testStreamSource_SOURCES = testStreamSource.cc
testMerge_SOURCES = testMerge.cc

bin_PROGRAMS = obs-print obs-info obs-skmeans \
               obs-concat obs-diff obs-stats obs-cat
#check_PROGRAMS = gmtest

#gmtest_SOURCES = gmtest.cc

obs_cat_SOURCES = obs-cat.cc
obs_stats_SOURCES = obs-stats.cc

obs_print_SOURCES = \
GMTK_ObsPrint.h GMTK_ObsPrint.cc \
GMTK_ObsKLT.h GMTK_ObsKLT.cc \
GMTK_ObsNorm.h GMTK_ObsNorm.cc \
GMTK_ObsGaussianNorm.h GMTK_ObsGaussianNorm.cc \
GMTK_ObsAddSil.h GMTK_ObsAddSil.cc

obs_info_SOURCES = GMTK_ObsInfo.cc

obs_skmeans_SOURCES = \
GMTK_ObsSkmeans.cc \
GMTK_ObsInitmg.h GMTK_ObsInitmg.cc \
GMTK_Kmeans.h GMTK_Kmeans.cc

obs_concat_SOURCES = GMTK_ObsConcat.cc

obs_diff_SOURCES = GMTK_ObsDiff.cc

