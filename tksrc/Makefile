#
#
# Makefile 
#   GMTK Makefile for GMTK's main toolkit source (tksrc) code.
# Do not run make directly using this directory. If you wish to compile only this module,
# from the current directory run the following command:
#
#    make -f ../Makefile MODULES=tksrc
#
# instead.
#
#
#
# $Header$
#
# 

SHELL = /bin/sh

# .EXPORT:
# .EXPORT: Linux AND !crane
# .EXPORT:


.SUFFIXES:
.SUFFIXES: .c .o .dis .cc

#
# Where sources are located.
#
srcdir = .
top_srcdir = /n/pato/home/bilmes/gmtoolkit/gmtk
VPATH = .

#
# default options
CC = gcc
CXX = g++
LD = $(CXX)
EXCCFLAGS=
EXCXXFLAGS=
EXLDFLAGS =  
OPTFLAGS =-g -O3
ANSI=-ansi
PEDANTIC=-pedantic  -Wno-long-long
WALL=-Wall
CCFLAGS = -g $(OPTFLAGS)  $(EXCCFLAGS) $(WALL) $(ANSI) $(PEDANTIC)
CXXFLAGS = -g $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

# extra optimization options for certain GMTK sources, used below.
XOPTFLAGS=-funroll-loops -fargument-noalias-global

AR = ar
LEX = flex
RANLIB = :
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 555
INSTALL_EXE = $(INSTALL) -m 755


prefix = .
exec_prefix = .
install_prefix = .

##  libdirs = -L$(prefix) -L$(prefix)/../miscSupport  -L$(prefix)/../pfile -L$(prefix)/../IEEEFloatingpoint -L$(prefix)/../featureFileIO 
libdirs = -L$(prefix) -L$(prefix)/../miscSupport  -L$(prefix)/../IEEEFloatingpoint -L$(prefix)/../featureFileIO
##  includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport -I$(prefix)/../pfile -I$(prefix)/../IEEEFloatingpoint -I$(prefix)/../featureFileIO 
includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport -I$(prefix)/../IEEEFloatingpoint -I$(prefix)/../featureFileIO
docdir = $(prefix)/man/man1

#
# Libraries.
#
# define only for sun
# SUNLIBS=-lsunmath
SUNLIBS=
LIBS =    -lGMTK -lgmtkio -lmiscSupport -lIEEEsupport -lfl $(SUNLIBS) -lm 
# LIBS =  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport -lfl $(SUNLIBS) -lm 
# LIBS =  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport $(SUNLIBS) -lm 
INCLUDES =  $(includedirs)
LDFLAGS= $(libdirs)
# Define PURE as purify to create purify versions.
# This can be done easily by doing 'make -DPURE=purify' 
PURE =
compile.cc = $(CXX) $(INCLUDES) $(CXXFLAGS)
compile.c = $(CC) $(INCLUDES) $(CCFLAGS)
link = $(PURE) $(LD) $(LDFLAGS) $(EXLDFLAGS) $(CXXFLAGS)
archive = $(AR) -r
compile.cc.link = $(PURE) $(LD)  $(INCLUDES) $(LDFLAGS) $(EXLDFLAGS) $(CXXFLAGS)

#----------------------------------------------------------------------------
# Default build rules.
#----------------------------------------------------------------------------
.c.o:
	-@rm -f $@
	$(compile.c) -c $<

.cc.o:
	-@rm -f $@
	$(compile.cc) -c $<

.cc.s:
	-@rm -f $@
	$(compile.cc) -S -c $<

#----------------------------------------------------------------------------
# Default target(s)
#----------------------------------------------------------------------------

all_exs = gmtkTie gmtkJT gmtkTime gmtkEMtrainNew gmtkKernel gmtkViterbiNew gmtkViterbi gmtkTriangulate gmtkParmConvert gmtkTFmerge gmtkDTindex gmtkNGramIndex gmtkTie #gmtkViz

all: ${all_exs}

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------

objs = GMTK_Dense1DPMF.o GMTK_Sparse1DPMF.o GMTK_MeanVector.o GMTK_DiagCovarVector.o GMTK_RealMatrix.o GMTK_DlinkMatrix.o GMTK_DirichletTable.o GMTK_DirichletPrior.o  GMTK_MDCPT.o GMTK_MSCPT.o GMTK_RealArray.o GMTK_CPT.o GMTK_MixtureCommon.o GMTK_Component.o GMTK_GaussianComponent.o GMTK_GMParms.o GMTK_FileParser.o lex.yy.o GMTK_MTCPT.o GMTK_Mixture.o GMTK_DiagGaussian.o GMTK_Dlinks.o GMTK_EMable.o GMTK_LinMeanCondDiagGaussian.o GMTK_LinMeanCondDiagGaussianOpt.o  GMTK_RngDecisionTree.o GMTK_NameCollection.o  GMTK_GMTemplate.o GMTK_JunctionTree.o GMTK_GraphicalModel.o GMTK_MaxClique.o GMTK_BoundaryTriangulate.o GMTK_Timer.o GMTK_Signals.o GMTK_PackCliqueValue.o GMTK_PackCliqueValueOpt.o GMTK_Vocab.o GMTK_NGramCPT.o GMTK_FNGramCPT.o GMTK_RV.o GMTK_ContRV.o GMTK_DiscRV.o GMTK_HidDiscRV.o GMTK_ObsDiscRV.o GMTK_ObsContRV.o  GMTK_ScPnShRV.o GMTK_SwDiscRV.o GMTK_SwRV.o GMTK_Sw_HidDiscRV.o GMTK_Sw_ObsContRV.o GMTK_Sw_ObsDiscRV.o GMTK_VECPT.o GMTK_NetworkFlow.o GMTK_NetworkFlowTriangulate.o GMTK_DiagGaussianOpt.o GMTK_LatticeADT.o GMTK_LatticeNodeCPT.o GMTK_LatticeEdgeCPT.o GMTK_Tie.o tieSupport.o GMTK_RVInfo.o GMTK_GammaComponent.o GMTK_BetaComponent.o GMTK_Partition.o GMTK_JT_Partition.o GMTK_JT_Partition.o GMTK_PartitionTables.o GMTK_PartitionStructures.o GMTK_JunctionTreeMessagePassing.o GMTK_JunctionTreeIslandMessagePassing.o






#----------------------------------------------------------------------------
# Other files used in the program. Keep this up to date.
#----------------------------------------------------------------------------


all_srcs =  GMTK_CPT.cc GMTK_Dlinks.cc GMTK_Dense1DPMF.cc GMTK_DiagCovarVector.cc GMTK_DiagGaussian.cc GMTK_DlinkMatrix.cc GMTK_GMParms.cc GMTK_MDCPT.cc GMTK_MSCPT.cc GMTK_MeanVector.cc GMTK_RealArray.cc GMTK_RealMatrix.cc GMTK_RngDecisionTree.cc GMTK_Sparse1DPMF.cc GMTK_DirichletTable.cc GMTK_DirichletPrior.cc  GMTK_FileParser.cc GMTK_FileTokenizer.lex GMTK_MTCPT.cc GMTK_MixtureCommon.cc GMTK_Component.cc GMTK_GaussianComponent.cc GMTK_Mixture.cc  gmtkEMtrainNew.cc gmtkKernel.cc GMTK_EMable.cc  gmtkViterbiNew.cc gmtkViterbi.cc GMTK_LinMeanCondDiagGaussian.cc GMTK_LinMeanCondDiagGaussianOpt.cc GMTK_NameCollection.cc GMTK_GraphicalModel.cc gmtkParmConvert.cc GMTK_GMTemplate.cc GMTK_JunctionTree.cc gmtkTriangulate.cc gmtkTFmerge.cc gmtkJT.cc gmtkTime.cc GMTK_MaxClique.cc GMTK_BoundaryTriangulate.cc GMTK_Timer.cc GMTK_Signals.cc gmtkDTindex.cc GMTK_PackCliqueValue.cc GMTK_PackCliqueValueOpt.cc GMTK_Vocab.cc GMTK_NGramCPT.cc GMTK_FNGramCPT.cc gmtkNGramIndex.cc GMTK_RV.cc GMTK_ContRV.cc GMTK_DiscRV.cc GMTK_HidDiscRV.cc GMTK_ObsContRV.cc GMTK_ObsDiscRV.cc GMTK_ScPnShRV.cc GMTK_ScPnSh_HidDiscRV.cc GMTK_SwDiscRV.cc GMTK_SwRV.cc GMTK_Sw_HidDiscRV.cc GMTK_Sw_ObsContRV.cc GMTK_Sw_ObsDiscRV.cc GMTK_VECPT.cc  GMTK_NetworkFlow.cc GMTK_NetworkFlowTriangulate.cc GMTK_DiagGaussianOpt.cc GMTK_LatticeADT.cc GMTK_LatticeNodeCPT.cc GMTK_LatticeEdgeCPT.cc gmtkTie.cc GMTK_Tie.cc tieSupport.cc GMTK_GammaComponent.cc GMTK_BetaComponent.cc GMTK_Partition.cc GMTK_JT_Partition.cc GMTK_PartitionTables.cc GMTK_PartitionStructures.cc GMTK_JunctionTreeMessagePassing.cc GMTK_JunctionTreeIslandMessagePassing.cc GMTK_MissingFeatureScaledDiagGaussian.cc GMTK_RVInfo.cc GMTK_WeightMatrix.cc
# gmtkUnroll.cc GMTK_RVInfo.cc gmtkViz.cc 



all_hdrs = GMTK_CPT.h GMTK_Constants.h GMTK_Dlinks.h GMTK_Dense1DPMF.h GMTK_DiagCovarVector.h GMTK_DiagGaussian.h GMTK_DiscreteRandomVariable.h GMTK_ContinuousRandomVariable.h GMTK_DlinkMatrix.h GMTK_EMable.h  GMTK_GMParms.h GMTK_MDCPT.h GMTK_MSCPT.h GMTK_MeanVector.h GMTK_RandomVariable.h GMTK_RealArray.h GMTK_RealMatrix.h GMTK_RngDecisionTree.h GMTK_Sparse1DPMF.h GMTK_DirichletTable.h GMTK_DirichletPrior.h  GMTK_FileParser.h GMTK_NamedObject.h GMTK_MTCPT.h GMTK_MixtureCommon.h GMTK_Component.h GMTK_GaussianComponent.h GMTK_Mixture.h GMTK_LinMeanCondDiagGaussian.h GMTK_ProgramDefaultParms.h  GMTK_NameCollection.h GMTK_GMTemplate.h GMTK_JunctionTree.h GMTK_GraphicalModel.h GMTK_ZeroScoreMixture.h GMTK_UnityScoreMixture.h GMTK_USCPT.h GMTK_RVInfo.h GMTK_FactorInfo.h GMTK_BoundaryTriangulate.h GMTK_MaxClique.h GMTK_Timer.h GMTK_Signals.h version.h GMTK_PackCliqueValue.h GMTK_SpaceManager.h GMTK_Vocab.h GMTK_NGramCPT.h GMTK_FNGramCPT.h GMTK_ContRV.h GMTK_DiscRV.h GMTK_HidDiscRV.h GMTK_ObsContRV.h GMTK_ObsDiscRV.h GMTK_RV.h GMTK_ScPnShRV.h GMTK_ScPnSh_HidDiscRV.h GMTK_ScPnSh_ObsContRV.h GMTK_ScPnSh_ObsDiscRV.h GMTK_ScPnSh_Sw_HidDiscRV.h GMTK_ScPnSh_Sw_ObsDiscRV.h GMTK_SimpleDiscreteRV.h GMTK_SwContRV.h GMTK_SwDiscRV.h GMTK_SwRV.h GMTK_Sw_HidDiscRV.h GMTK_Sw_ObsContRV.h GMTK_Sw_ObsDiscRV.h GMTK_VECPT.h GMTK_NetworkFlow.h GMTK_NetworkFlowTriangulate.h GMTK_LatticeADT.h GMTK_LatticeNodeCPT.h GMTK_LatticeEdgeCPT.h GMTK_Tie.h tieSupport.h GMTK_GammaComponent.h GMTK_BetaComponent.h GMTK_Partition.h GMTK_JT_Partition.h GMTK_PartitionTables.h GMTK_PartitionStructures.h GMTK_Arguments.h GMTK_DiscRVType.h GMTK_MissingFeatureScaledDiagGaussian.h GMTK_ScPnSh_Sw_ObsContRV.h GMTK_WeightMatrix.h





all_objs = $(objs)

#----------------------------------------------------------------------------
# Special case rules for .cc -> .o
#----------------------------------------------------------------------------

GMTK_DiagGaussianOpt.o: GMTK_DiagGaussianOpt.cc
	-@rm -f $@
	$(compile.cc) $(XOPTFLAGS) -c $<

GMTK_LinMeanCondDiagGaussianOpt.o: GMTK_LinMeanCondDiagGaussianOpt.cc
	-@rm -f $@
	$(compile.cc) $(XOPTFLAGS) -c $<


GMTK_PackCliqueValueOpt.o: GMTK_PackCliqueValueOpt.cc
	-@rm -f $@
	$(compile.cc) $(XOPTFLAGS) -c $<

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

.NOEXPORT:

programname: $(objs)
	$(link) -Wall -g  -o $@ $(bmm_objs) $@.o $(LIBS)


libGMTK.a: $(objs)
	rm -f $@
	$(archive) $@ $(objs)


gmtkEMtrainNew: gmtkEMtrainNew.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkEMtrainNew.o $(LIBS)

gmtkKernel: gmtkKernel.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkKernel.o $(LIBS)


#gmtkUnroll: gmtkUnroll.o  libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
#	-@rm -f $@
#	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkUnroll.o $(LIBS)

gmtkJT: gmtkJT.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkJT.o $(LIBS)

gmtkTime: gmtkTime.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTime.o $(LIBS)


gmtkTriangulate: gmtkTriangulate.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTriangulate.o $(LIBS)


gmtkTFmerge: gmtkTFmerge.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTFmerge.o $(LIBS)


gmtkViterbiNew: gmtkViterbiNew.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkViterbiNew.o $(LIBS)

gmtkViterbi: gmtkViterbi.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkViterbi.o $(LIBS)


gmtkViterbiNew3: gmtkViterbiNew3.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkViterbiNew3.o $(LIBS)

gmtkParmConvert: gmtkParmConvert.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkParmConvert.o $(LIBS)

gmtkDTindex: gmtkDTindex.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkDTindex.o $(LIBS)

gmtkNGramIndex: gmtkNGramIndex.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkNGramIndex.o $(LIBS)

gmtkTie: gmtkTie.o libGMTK.a ../miscSupport/libmiscSupport.a .NOEXPORT
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTie.o $(LIBS)

gmtkViz.o: gmtkViz.cc 
	-@rm -f $@
	$(CXX) $(INCLUDES) -g $(OPTFLAGS)  $(EXCFLAGS) -Wall   -Wno-deprecated `wx-config --cxxflags` -Wall -g -c gmtkViz.cc

gmtkViz: gmtkViz.o libGMTK.a
	-@rm -f $@
	$(compile.cc.link) -Wall -g -o $@ gmtkViz.o $(LIBS) `wx-config --libs` `wx-config --cxxflags` 

#----------------------------------------------------------------------------
# Test driver programs (defines a -DMAIN to test each .cc file)
#----------------------------------------------------------------------------

testMDCPT: GMTK_MDCPT.cc GMTK_MDCPT.h  GMTK_CPT.o
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_MDCPT.cc GMTK_CPT.o $(LIBS)


testMTCPT: GMTK_MTCPT.cc GMTK_MTCPT.h  $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_MTCPT.cc $(objs) $(LIBS)


testHash: GMTK_Hash.cc GMTK_Hash.h 
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_Hash.cc $(LIBS)

testHashTree: GMTK_HashTree.cc GMTK_HashTree.h 
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_HashTree.cc $(LIBS)

testRngDecisionTree: GMTK_RngDecisionTree.cc GMTK_RngDecisionTree.h GMTK_RV.cc GMTK_RV.h GMTK_DiscRV.cc GMTK_DiscRV.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_RngDecisionTree.cc GMTK_RV.o GMTK_DiscRV.o $(LIBS)

testPackCliqueValue: GMTK_PackCliqueValue.cc  GMTK_PackCliqueValue.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_PackCliqueValue.cc $(LIBS)

testGMParms: GMTK_GMParms.cc GMTK_GMParms.h $(objs)
	echo $(objs) > /tmp/gmparms_depts
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_GMParms.cc `sed s,GMTK_GMParms.o,, < /tmp/gmparms_depts` $(LIBS)


lex.yy: lex.yy.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ lex.yy.cc -lfl

lex.yy.cc: GMTK_FileTokenizer.lex
	${LEX} -B GMTK_FileTokenizer.lex
	mv lex.yy.c lex.yy.cc	

testFileParser: GMTK_FileParser.cc GMTK_FileParser.h lex.yy.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_FileParser.cc lex.yy.o $(objs) $(LIBS) -lfl


#----------------------------------------------------------------------------
# Installation
#----------------------------------------------------------------------------

install-docs:

uninstall-docs:

install: $(all_exs)
	(for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i; done)

installd: $(all_exs)
	(STR=`date +%a_%b_%d_%Y_%k:%M | sed -e 's, ,,g'`; echo Installing with postfix $$STR; for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i.$$STR; strip $(install_prefix)/bin/$$i.$$STR; done)

installbackexs:
	(STR=`date +%a_%b_%d_%Y_%k:%M | sed -e 's, ,,g'`; echo Installing with postfix $$STR; for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i.$$STR; strip $(install_prefix)/bin/$$i.$$STR; done)

uninstall:
	rm -f $(exec_prefix)/bin/foo

#----------------------------------------------------------------------------
# TODO:Configure targets.
#----------------------------------------------------------------------------
#
#$(srcdir)/configure: $(srcdir)/configure.in $(srcdir)/aclocal.m4
#	cd $(srcdir) && autoconf
#
#Makefile: $(srcdir)/Makefile.in config.status
#	./config.status
#     
#config.status: $(srcdir)/configure
#	./config.status --recheck
#
#config_files = $(srcdir)/configure.in $(srcdir)/aclocal.m4 \
#		$(srcdir)/Makefile.in
# Include this for now until we get configure running.
config_files = Makefile TODO README notes
#
#----------------------------------------------------------------------------
# TODO: Documentation.
#----------------------------------------------------------------------------
#
#doc_files = $(srcdir)/INSTALL $(srcdir)/NEWS $(srcdir)/TODO \
#		$(srcdir)/README
#
#----------------------------------------------------------------------------
# RCS management.
#----------------------------------------------------------------------------

# VERSION is symbolic name of version for RCS targets.
# CHANGES is message to use when making bulk check-ins.

rcs_files = $(all_srcs) $(all_hdrs) $(doc_files) $(config_files)

ci:
	ci -u -t- "-m$(CHANGES)" $(rcs_files)

co:
	co -l $(rcs_files)

cio:
	ci -l -t- "-m$(CHANGES)" $(rcs_files)

# Check out latest copy.
col:
	co $(rcs_files)

# Check out a previous version.
cov:
	co -r$(VERSION) $(rcs_files)

rcsdiff:
	rcsdiff $(rcs_files)

# Make sure all files are checked in before creating new version.
ensure-ci:
	@if [ `rlog -L -l -R $(rcs_files) | wc -w` -ne 0 ]; \
	then \
		echo "ERROR: The following files are not checked in:"; \
		rlog -L -l -h $(rcs_files) | egrep 'Working file|bilmes' ; \
		exit 1; \
	fi

#
#----------------------------------------------------------------------------
# CVS management.
#----------------------------------------------------------------------------

TODO-dev:
	cvs ci -m'update' TODO; cvs tag -F development TODO

update:
	cvs update

# make current version a development version
dev:
	cvs tag -F development

diff:
	cvs diff

# to make a new version, do something like:
#	make VERSION='June12Working' version	
version :
	cvs tag -F "$(VERSION)"



#----------------------------------------------------------------------------
# Generate TAGS file for emacs editing.
#----------------------------------------------------------------------------

TAGS: $(all_srcs) $(all_hdrs)
	etags --output=$(srcdir)/TAGS $(all_srcs) $(all_hdrs)

#----------------------------------------------------------------------------
# Count lines of code and number of C statements.
#----------------------------------------------------------------------------

counts:
	wc $(all_srcs) $(all_hdrs)
	cat $(all_srcs) $(all_hdrs) | tr -cd ';' | wc -c

#----------------------------------------------------------------------------
# Garbage disposal.
#----------------------------------------------------------------------------

target_junk = $(all_objs) $(all_exs) $(lib_name)

clean:
	rm -rf $(target_junk) *.o

exsclean: 
	rm -rf $(all_exs)

# distclean: clean
# 	rm -f config.log config.status config.cache Makefile

maintainer_junk = $(srcdir)/TAGS

maintainer-clean: distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."
	rm -rf $(maintainer_junk)

depend: 
	${CXX} ${INCLUDES} -MM ${CXXFLAGS} ${all_srcs} > depends.make

depends: depend

##################
####  DEPENDS ####
include depends.make
