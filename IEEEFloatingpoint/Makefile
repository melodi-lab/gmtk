#
# Makefile 
#   for misc support
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
top_srcdir = .
VPATH = .

#
# Compilng programs to use.
#
CXX = g++
CC = gcc
LD = $(CXX)
AR = ar
RANLIB = :
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 555
INSTALL_EXE = $(INSTALL) -m 555

#
# Program options.
#
# OPTFLAGS = -O
OPTFLAGS = -O -DHAVE_NONSTANDARD_ARITHMETIC

# extra flags via command line (e.g., -pg).
EXCFLAGS =  
EXLDFLAGS =  
CXXFLAGS = -g $(OPTFLAGS) $(EXCFLAGS) -Wall -ansi -pedantic
CCFLAGS = -g $(OPTFLAGS)  $(EXCFLAGS) -Wall -ansi -pedantic
LDFLAGS = 
DEFS = 


prefix = .
exec_prefix = .
install_prefix = .

libdirs = -L$(prefix)
includedirs = -I$(srcdir) -I$(srcdir)/../miscSupport
docdir = $(prefix)/man/man1

#
# Libraries.
#
LIBS = -lm
INCLUDES =  $(includedirs)
LDFLAGS= $(libdirs)
compile.cc = $(CXX) $(DEFS) $(INCLUDES) $(CCFLAGS)
compile.c = $(CC) $(DEFS) $(INCLUDES) $(CXXFLAGS)
pure =
link = $(pure) $(LD) $(LDFLAGS) $(EXLDFLAGS) $(CXXFLAGS)
archive = $(AR) -r
compile.cc.link = $(PURE) $(LD)  $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

#----------------------------------------------------------------------------
# Default build rules.
#----------------------------------------------------------------------------
.c.o:
	$(compile.c) -c $<

.cc.o:
	$(compile.cc) -c $<

.cc.s:
	$(compile.cc) -S -c $<

#----------------------------------------------------------------------------
# Default target(s)
#----------------------------------------------------------------------------

all: libIEEEsupport.a

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------

objs = ieeeFPsetup.o

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

libIEEEsupport.a: $(objs)
	rm -f $@
	$(archive) $@ $(objs)

clean_depends:
	rm depends.make
	touch depends.make
	make depends

solaris:
	rm -f ieeeFPsetup.h
	rm -f ieeeFPsetup.cc
	ln -s ieeeFPsetup_solaris.h ieeeFPsetup.h
	ln -s ieeeFPsetup_solaris.cc ieeeFPsetup.cc


linux:
	rm -f ieeeFPsetup.h
	rm -f ieeeFPsetup.cc
	ln -s ieeeFPsetup_linux.h ieeeFPsetup.h
	ln -s ieeeFPsetup_linux.cc ieeeFPsetup.cc


ibm:
	rm -f ieeeFPsetup.h
	rm -f ieeeFPsetup.cc
	ln -s ieeeFPsetup_ibm.h ieeeFPsetup.h
	ln -s ieeeFPsetup_ibm.cc ieeeFPsetup.cc

cygwin:
	rm -f ieeeFPsetup.h
	rm -f ieeeFPsetup.cc
	ln -s ieeeFPsetup_cygwin.h ieeeFPsetup.h
	ln -s ieeeFPsetup_cygwin.cc ieeeFPsetup.cc


#----------------------------------------------------------------------------
# Test programs for each .cc file
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Gather stuff together.
#----------------------------------------------------------------------------

all_srcs = ieeeFPsetup.cc

all_hdrs = ieeeFPsetup.h

all_objs = 

all_exs = 

#----------------------------------------------------------------------------
# installation
#----------------------------------------------------------------------------

install-docs:

uninstall-docs:

install: $(all_exs)
	(for i in $(all_exs); do $(INSTALL_EXE) $$i $(install_prefix)/bin/$$i; done)

uninstall:
	rm -f $(exec_prefix)/bin/pfile_stats

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
		rlog -L -l -R $(rcs_files); \
		exit 1; \
	fi

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


##################
####  DEPENDS ####
# To make dependencies, first do:
#     make depend
# if the file doesn't exist, do a 'touch depends.make' first.
depends: depend

# include the dependencies
include depends.make


