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
INSTALL = /usr/ucb/install -c
INSTALL_DATA = ${INSTALL} -m 555
INSTALL_EXE = $(INSTALL) -m 555

#
# Program options.
#
# OPTFLAGS = -O2 -DNDEBUG
OPTFLAGS=-O3

EXFLAGS =  # extra flags via command line (e.g., -pg).
CXXFLAGS = -g $(OPTFLAGS) $(EXFLAGS) -Wall -ansi -pedantic
CCFLAGS = -g $(OPTFLAGS)  $(EXFLAGS) -Wall -ansi -pedantic
LDFLAGS = 
DEFS = 


prefix = .
exec_prefix = .
install_prefix = .

libdirs = -L$(prefix) -L$(prefix)/../miscSupport -L$(prefix)/../lowLevelFileIO -L$(prefix)/../pfile -L$(prefix)/../IEEEFloatingpoint -L$(prefix)/../featureFileIO -L/g/rcs/sw/lib
includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport -I$(prefix)/../lowLevelFileIO -I$(prefix)/../pfile -I$(prefix)/../IEEEFloatingpoint -I$(prefix)/../featureFileIO 
docdir = $(prefix)/man/man1

#
# Libraries.
#
SUNLIBS=-lsunmath
LIBS = -lfileParser  -lmiscSupport  -lgmtkio -lpfile -lIEEEsupport -lfl $(SUNLIBS) -lm 
INCLUDES =  $(includedirs)
LDFLAGS= $(libdirs)
# Define PURE as purify to create purify versions.
# This can be done easily by doing 'make -DPURE=purify' 
PURE =
compile.cc = $(CXX) $(DEFS) $(INCLUDES) $(CCFLAGS)
compile.c = $(CC) $(DEFS) $(INCLUDES) $(CXXFLAGS)
link = $(PURE) $(LD) $(LDFLAGS) $(CXXFLAGS)
archive = $(AR) -r
compile.cc.link = $(PURE) $(LD)  $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

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

all_exs = gmtkEMtrain gmtkViterbi gmtkM2T gmtkB2A gmtkScore

all: ${all_exs}

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------

objs = GMTK_Dense1DPMF.o GMTK_Sparse1DPMF.o GMTK_MeanVector.o GMTK_DiagCovarVector.o GMTK_RealMatrix.o GMTK_DlinkMatrix.o GMTK_WeightMatrix.o GMTK_MDCPT.o GMTK_MSCPT.o GMTK_RealArray.o GMTK_CPT.o GMTK_MixGaussiansCommon.o GMTK_GaussianComponent.o GMTK_GMParms.o GMTK_RandomVariable.o GMTK_DiscreteRandomVariable.o  GMTK_SimpleDiscreteRV.o GMTK_ContinuousRandomVariable.o GMTK_FileParser.o lex.yy.o GMTK_MTCPT.o GMTK_GM.o  GMTK_CliqueChain.o GMTK_Clique.o GMTK_MixGaussians.o GMTK_DiagGaussian.o GMTK_Dlinks.o GMTK_EMable.o GMTK_LinMeanCondDiagGaussian.o GMTK_Hash.o

#----------------------------------------------------------------------------
# Other files used in the program. Keep this up to date.
#----------------------------------------------------------------------------

all_srcs = gmtkScore.cc BN1.cc GMTK_CPT.cc GMTK_Clique.cc GMTK_CliqueChain.cc GMTK_Dlinks.cc GMTK_Dense1DPMF.cc GMTK_DiagCovarVector.cc GMTK_DiagGaussian.cc GMTK_DiscreteRandomVariable.cc GMTK_ContinuousRandomVariable.cc GMTK_DlinkMatrix.cc GMTK_GM.cc GMTK_GMParms.cc GMTK_MDCPT.cc GMTK_MSCPT.cc GMTK_MeanVector.cc GMTK_RandomVariable.cc GMTK_RealArray.cc GMTK_RealMatrix.cc GMTK_RngDecisionTree.cc GMTK_Sparse1DPMF.cc GMTK_WeightMatrix.cc GMTK_FileParser.cc GMTK_FileTokenizer.lex GMTK_MTCPT.cc GMTK_MixGaussiansCommon.cc GMTK_GaussianComponent.cc GMTK_MixGaussians.cc gmtkEMtrain.cc GMTK_EMable.cc gmtkViterbi.cc GMTK_LinMeanCondDiagGaussian.cc GMTK_Hash.cc

all_hdrs = GMTK_CPT.h GMTK_Clique.h GMTK_CliqueChain.h GMTK_Constants.h GMTK_Dlinks.h GMTK_Dense1DPMF.h GMTK_DiagCovarVector.h GMTK_DiagGaussian.h GMTK_DiscreteRandomVariable.h GMTK_ContinuousRandomVariable.h GMTK_DlinkMatrix.h GMTK_EMable.h GMTK_GM.h GMTK_GMParms.h GMTK_MDCPT.h GMTK_MSCPT.h GMTK_MeanVector.h GMTK_RandomVariable.h GMTK_RealArray.h GMTK_RealMatrix.h GMTK_RngDecisionTree.h GMTK_Sparse1DPMF.h GMTK_WeightMatrix.h GMTK_FileParser.h GMTK_NamedObject.h GMTK_MTCPT.h GMTK_MixGaussiansCommon.h GMTK_GaussianComponent.h GMTK_MixGaussians.h GMTK_LinMeanCondDiagGaussian.h GMTK_Hash.h

all_objs = $(objs)

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

programname: $(objs)
	$(link) -Wall -g  -o $@ $(bmm_objs) $@.o $(LIBS)

gmtkEMtrain: gmtkEMtrain.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkEMtrain.o $(objs) $(LIBS)

gmtkViterbi: gmtkViterbi.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkViterbi.o $(objs) $(LIBS)

gmtkScore: gmtkScore.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkScore.o $(objs) $(LIBS)

gmtkM2T: gmtkM2T.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkM2T.o $(objs) $(LIBS)


gmtkB2A: gmtkB2A.o $(objs)
	-@rm -f $@
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ gmtkB2A.o $(objs) $(LIBS)

#----------------------------------------------------------------------------
# Test driver programs (defines a -DMAIN to test each .cc file)
#----------------------------------------------------------------------------

testMDCPT: GMTK_MDCPT.cc GMTK_MDCPT.h  GMTK_CPT.o
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_MDCPT.cc GMTK_CPT.o $(LIBS)


testMTCPT: GMTK_MTCPT.cc GMTK_MTCPT.h  $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_MTCPT.cc $(objs) $(LIBS)


testHash: GMTK_Hash.cc GMTK_Hash.h 
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_Hash.cc $(LIBS)


testRngDecisionTree: GMTK_RngDecisionTree.cc GMTK_RngDecisionTree.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_RngDecisionTree.cc $(LIBS)

testGMParms: GMTK_GMParms.cc GMTK_GMParms.h $(objs)
	echo $(objs) > /tmp/gmparms_depts
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ GMTK_GMParms.cc `sed s,GMTK_GMParms.o,, < /tmp/gmparms_depts` $(LIBS)


testBN1: testBN1.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ testBN1.o $(objs) $(LIBS)

BN1: BN1.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ BN1.o $(objs) $(LIBS)

BNU: BNU.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ BNU.o $(objs) $(LIBS)

testBN2: testBN2.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ testBN2.o $(objs) $(LIBS)

BN2: BN2.cc BN2.o GMTK_GM.cc GMTK_GM.o GMTK_Clique.cc GMTK_Clique.o \
GMTK_CliqueChain.cc GMTK_CliqueChain.o GMTK_RandomVariable.cc \
GMTK_RandomVariable.o GMTK_DiscreteRandomVariable.cc \
GMTK_DiscreteRandomVariable.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ BN2.cc GMTK_GM.cc \
GMTK_Clique.cc GMTK_CliqueChain.cc GMTK_RandomVariable.cc \
GMTK_DiscreteRandomVariable.cc $(objs) $(LIBS)

lex.yy: lex.yy.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ lex.yy.cc -lfl

lex.yy.cc: GMTK_FileTokenizer.lex
	${LEX} GMTK_FileTokenizer.lex
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
config_files = Makefile TODO README data_items_needed
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
	${CXX} ${INCLUDES} -MM ${CXXFLAGS} *.cc *.c > depends.make

depends: depend

##################
####  DEPENDS ####
include depends.make


