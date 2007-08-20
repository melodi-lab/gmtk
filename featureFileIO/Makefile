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
CC =  gcc
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

# extra flags via command line (e.g., -pg).
EXCFLAGS =  
EXLDFLAGS =  
CXXFLAGS = -g $(OPTFLAGS) $(EXCFLAGS) -Wall -ansi -DMAIN
CCFLAGS = -g $(OPTFLAGS)  $(EXCFLAGS) -Wall -ansi -DMAIN
LDFLAGS = 
DEFS = 

# doesn't appear to work with pmake, temporarly remove. JB:12/29/04.
# i f d e f DEBUG
# CXXFLAGS := $(CXXFLAGS) -DDEBUG='$(DEBUG)'
# e n d i f


prefix = .
exec_prefix = .
install_prefix = .

libdirs = -L$(prefix) -L$(prefix)/../miscSupport -L$(prefix)/../IEEEFloatingpoint -L/g/rcs/sw/lib
includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport -I$(prefix)/../IEEEFloatingpoint
docdir = $(prefix)/man/man1


LIBS = -lfileParser  -lmiscSupport  -lIEEEsupport -lfl -lsunmath -lm 
GMTEST_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm 
OBS_PRINT_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm
OBS_INFO_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm
OBS_SKMEANS_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm
OBS_CONCAT_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm
OBS_DIFF_LIBS = -lmiscSupport  -lIEEEsupport -lfl -lm


INCLUDES =  $(includedirs)
LDFLAGS= $(libdirs)
# Define PURE as purify to create purify versions.
# This can be done easily by doing 'make -DPURE=purify' 
PURE =
compile.cc = $(CXX) $(DEFS) $(INCLUDES) $(CXXFLAGS)
compile.c = $(CC) $(DEFS) $(INCLUDES) $(CCFLAGS)
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

all_exs = libgmtkio.a obs-print obs-info obs-skmeans obs-concat obs-diff

all: ${all_exs}

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------

objs = GMTK_ObservationMatrix.o GMTK_Stream.o GMTK_DataTransformations.o GMTK_WordOrganization.o GMTK_ObsStats.o pfile.o vbyteswapping.o

#----------------------------------------------------------------------------
# Other files used in the program. Keep this up to date.
#----------------------------------------------------------------------------

all_srcs = GMTK_ObservationMatrix.cc GMTK_Stream.cc GMTK_DataTransformations.cc GMTK_WordOrganization.cc  GMTK_ObsStats.cc vbyteswapping.cc

all_hdrs = GMTK_Stream.h GMTK_ObservationMatrix.h GMTK_WordOrganization.h  GMTK_ObsStats.h vbyteswapping.h

all_objs = $(objs)

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

gmtest: $(objs) gmtest.o
	rm -f $@
	$(link) -Wall -g  -o $@ $@.o $(objs) $(GMTEST_LIBS)


libgmtkio.a: $(objs)
	rm -f $@
	$(archive) $@ $(objs)

obs-print: $(objs) GMTK_ObsPrint.o GMTK_ObsKLT.o GMTK_ObsNorm.o GMTK_ObsGaussianNorm.o GMTK_ObsAddSil.o
	$(link) -Wall -o $@ GMTK_ObsPrint.o  GMTK_ObsKLT.o GMTK_ObsNorm.o GMTK_ObsGaussianNorm.o GMTK_ObsAddSil.o $(objs) $(OBS_PRINT_LIBS)

obs-info: $(objs) GMTK_ObsInfo.o
	$(link) -Wall -o $@ GMTK_ObsInfo.o $(objs) $(OBS_INFO_LIBS)

obs-skmeans: $(objs) GMTK_ObsSkmeans.o GMTK_ObsInitmg.o GMTK_Kmeans.o
	$(link) -Wall -o $@ GMTK_ObsSkmeans.o GMTK_ObsInitmg.o GMTK_Kmeans.o $(objs) $(OBS_SKMEANS_LIBS)

obs-concat: $(objs) GMTK_ObsConcat.o
	 $(link) -Wall -o $@ GMTK_ObsConcat.o $(objs) $(OBS_CONCAT_LIBS)

obs-diff: $(objs) GMTK_ObsDiff.o
	$(link) -Wall -o $@ GMTK_ObsDiff.o $(objs) $(OBS_DIFF_LIBS)

#----------------------------------------------------------------------------
# Test driver programs (defines a -DMAIN to test each .cc file)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Installation
#----------------------------------------------------------------------------

install-docs:

uninstall-docs:

install: $(all_exs)
	(for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i; done)

installd: $(all_exs)
	(for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i.`date +%a%b%d_%k_%Y`; done)

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
config_files = Makefile
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
