#
# Makefile 
#   General Makefile for GMTK
#   J. Bilmes
#
#
# $Header$
#
# 

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .o .dis .cc

#
# Where sources are located.
#
srcdir = .
top_srcdir = /n/pato/home/bilmes/gmtoolkit/gmtk
VPATH = .

#
# Compilng programs to use.
#
CXX = g++
CC = gcc
LD = $(CXX)
AR = ar
LEX = flex
RANLIB = :
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 555
INSTALL_EXE = $(INSTALL) -m 755

#
# Program options.
#
# OPTFLAGS = -O2 -DNDEBUG
OPTFLAGS=-O3
# OPTFLAGS=

# extra flags via command line (e.g., -pg).
EXCFLAGS =  
EXLDFLAGS =  
CXXFLAGS = -g $(OPTFLAGS) $(EXCFLAGS) -Wall  -Wno-deprecated -ansi -pedantic
CCFLAGS = -g $(OPTFLAGS)  $(EXCFLAGS) -Wall   -Wno-deprecated -ansi -pedantic
LDFLAGS = 
DEFS = 


prefix = .
exec_prefix = .
install_prefix = .

libdirs = -L$(prefix) -L$(prefix)/../miscSupport  -L$(prefix)/../pfile -L$(prefix)/../IEEEFloatingpoint -L$(prefix)/../featureFileIO -L/g/rcs/sw/lib
includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport -I$(prefix)/../pfile -I$(prefix)/../IEEEFloatingpoint -I$(prefix)/../featureFileIO 
docdir = $(prefix)/man/man1

#
# Libraries.
#
# define only for sun
# SUNLIBS=-lsunmath
SUNLIBS=
LIBS =  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport -lfl $(SUNLIBS) -lm 
# LIBS =  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport -lfl $(SUNLIBS) -lm 
# LIBS =  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport $(SUNLIBS) -lm 
INCLUDES =  $(includedirs)
LDFLAGS= $(libdirs)
# Define PURE as purify to create purify versions.
# This can be done easily by doing 'make -DPURE=purify' 
PURE =
compile.cc = $(CXX) $(DEFS) $(INCLUDES) $(CCFLAGS)
compile.c = $(CC) $(DEFS) $(INCLUDES) $(CXXFLAGS)
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

all_exs = gmtkEMtrain gmtkViterbi gmtkParmConvert gmtkScore gmtkSample gmtkTriangulate gmtkTFmerge gmtkJT gmtkUnroll gmtkDTindex

all: ${all_exs}

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------


objs = GMTK_Dense1DPMF.o GMTK_Sparse1DPMF.o GMTK_MeanVector.o GMTK_DiagCovarVector.o GMTK_RealMatrix.o GMTK_DlinkMatrix.o GMTK_WeightMatrix.o GMTK_MDCPT.o GMTK_MSCPT.o GMTK_RealArray.o GMTK_CPT.o GMTK_MixtureCommon.o GMTK_Component.o GMTK_GaussianComponent.o GMTK_GMParms.o GMTK_RandomVariable.o GMTK_DiscreteRandomVariable.o  GMTK_ContinuousRandomVariable.o GMTK_FileParser.o lex.yy.o GMTK_MTCPT.o GMTK_GM.o  GMTK_CliqueChain.o GMTK_Clique.o GMTK_Mixture.o GMTK_DiagGaussian.o GMTK_Dlinks.o GMTK_EMable.o GMTK_LinMeanCondDiagGaussian.o GMTK_Hash.o GMTK_RngDecisionTree.o GMTK_NameCollection.o  GMTK_GMTemplate.o GMTK_JunctionTree.o GMTK_GraphicalModel.o GMTK_MaxClique.o GMTK_BoundaryTriangulate.o GMTK_Timer.o GMTK_Signals.o GMTK_PackCliqueValue.o




#----------------------------------------------------------------------------
# Other files used in the program. Keep this up to date.
#----------------------------------------------------------------------------


all_srcs = gmtkSample.cc gmtkScore.cc gmtkUnroll.cc GMTK_CPT.cc GMTK_Clique.cc GMTK_CliqueChain.cc GMTK_Dlinks.cc GMTK_Dense1DPMF.cc GMTK_DiagCovarVector.cc GMTK_DiagGaussian.cc GMTK_DiscreteRandomVariable.cc GMTK_ContinuousRandomVariable.cc GMTK_DlinkMatrix.cc GMTK_GM.cc GMTK_GMParms.cc GMTK_MDCPT.cc GMTK_MSCPT.cc GMTK_MeanVector.cc GMTK_RandomVariable.cc GMTK_RealArray.cc GMTK_RealMatrix.cc GMTK_RngDecisionTree.cc GMTK_Sparse1DPMF.cc GMTK_WeightMatrix.cc GMTK_FileParser.cc GMTK_FileTokenizer.lex GMTK_MTCPT.cc GMTK_MixtureCommon.cc GMTK_Component.cc GMTK_GaussianComponent.cc GMTK_Mixture.cc gmtkEMtrain.cc GMTK_EMable.cc gmtkViterbi.cc GMTK_LinMeanCondDiagGaussian.cc GMTK_Hash.cc GMTK_NameCollection.cc GMTK_GraphicalModel.cc gmtkParmConvert.cc GMTK_GMTemplate.cc GMTK_JunctionTree.cc gmtkTriangulate.cc gmtkTFmerge.cc gmtkJT.cc GMTK_MaxClique.cc GMTK_BoundaryTriangulate.cc GMTK_Timer.cc GMTK_Signals.cc gmtkDTindex.cc GMTK_PackCliqueValue.cc

all_hdrs = GMTK_CPT.h GMTK_Clique.h GMTK_CliqueChain.h GMTK_Constants.h GMTK_Dlinks.h GMTK_Dense1DPMF.h GMTK_DiagCovarVector.h GMTK_DiagGaussian.h GMTK_DiscreteRandomVariable.h GMTK_ContinuousRandomVariable.h GMTK_DlinkMatrix.h GMTK_EMable.h GMTK_GM.h GMTK_GMParms.h GMTK_MDCPT.h GMTK_MSCPT.h GMTK_MeanVector.h GMTK_RandomVariable.h GMTK_RealArray.h GMTK_RealMatrix.h GMTK_RngDecisionTree.h GMTK_Sparse1DPMF.h GMTK_WeightMatrix.h GMTK_FileParser.h GMTK_NamedObject.h GMTK_MTCPT.h GMTK_MixtureCommon.h GMTK_Component.h GMTK_GaussianComponent.h GMTK_Mixture.h GMTK_LinMeanCondDiagGaussian.h GMTK_Hash.h GMTK_ProgramDefaultParms.h  GMTK_NameCollection.h GMTK_GMTemplate.h GMTK_JunctionTree.h GMTK_GraphicalModel.h GMTK_ZeroScoreMixture.h GMTK_UnityScoreMixture.h GMTK_USCPT.h GMTK_RVInfo.h GMTK_BoundaryTriangulate.h GMTK_MaxClique.h GMTK_Timer.h GMTK_Signals.h version.h GMTK_PackCliqueValue.h

all_objs = $(objs)

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

programname: $(objs)
	$(link) -Wall -g  -o $@ $(bmm_objs) $@.o $(LIBS)

gmtkEMtrain: gmtkEMtrain.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkEMtrain.o $(objs) $(LIBS)

gmtkUnroll: gmtkUnroll.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkUnroll.o $(objs) $(LIBS)

gmtkJT: gmtkJT.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkJT.o $(objs) $(LIBS)

gmtkTriangulate: gmtkTriangulate.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTriangulate.o $(objs) $(LIBS)


gmtkTFmerge: gmtkTFmerge.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkTFmerge.o $(objs) $(LIBS)


gmtkViterbi: gmtkViterbi.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkViterbi.o $(objs) $(LIBS)

gmtkScore: gmtkScore.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkScore.o $(objs) $(LIBS)

gmtkSample: gmtkSample.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkSample.o $(objs) $(LIBS)

gmtkParmConvert: gmtkParmConvert.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkParmConvert.o $(objs) $(LIBS)

gmtkDTindex: gmtkDTindex.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkDTindex.o $(objs) $(LIBS)

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

testRngDecisionTree: GMTK_RngDecisionTree.cc GMTK_RngDecisionTree.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_RngDecisionTree.cc $(LIBS)

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
	(STR=`date +%a%b%d_%k_%Y | sed -e 's, ,,g'`; echo Installing with postfix $$STR; for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i.$$STR; strip $(install_prefix)/bin/$$i.$$STR; done)

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
		rlog -L -l -h $(rcs_files) | egrep 'Working file|bilmes|zweig' ; \
		exit 1; \
	fi

# to make a new version, do something like:
#	make VERSION='June12Working' version	
version : ensure-ci
	rcs "-n$(VERSION):" $(rcs_files)

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

distclean: clean
	rm -f config.log config.status config.cache Makefile

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


