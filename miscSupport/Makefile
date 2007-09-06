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
OPTFLAGS = -O3
# OPTFLAGS =

# extra flags via command line (e.g., -pg).
EXCFLAGS =  
EXLDFLAGS =  
# CXXFLAGS = -g $(OPTFLAGS) $(EXCFLAGS) -Wall -ansi -pedantic
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
compile.cc = $(CXX) $(DEFS) $(INCLUDES) $(CXXFLAGS) $(HASH)
compile.c = $(CC) $(DEFS) $(INCLUDES) $(CCFLAGS) $(HASH)
pure =
link = $(pure) $(LD) $(LDFLAGS) $(EXLDFLAGS) $(CXXFLAGS)
archive = $(AR) -r
compile.cc.link = $(PURE) $(LD)  $(INCLUDES) $(LDFLAGS) $(EXLDFLAGS) $(CXXFLAGS) $(HASH)

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

all: libmiscSupport.a 

#----------------------------------------------------------------------------
# Define the standard library modules
#----------------------------------------------------------------------------

objs = error.o bp_range.o range.o arguments.o logp.o rand.o sArray.o cArray.o general.o lineqsolve.o debug.o fileParser.o hash_abstract.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o eig.o bin_search.o

#----------------------------------------------------------------------------
# Main build programs.
#----------------------------------------------------------------------------

libmiscSupport.a: $(objs)
	rm -f $@
	$(archive) $@ $(objs)


#----------------------------------------------------------------------------
# Test programs for each .cc file
#----------------------------------------------------------------------------

testError: error.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ error.cc $(objs)

testDebug: debug.cc debug.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ debug.cc $(objs)


testBPRange: bp_range.cc error.o
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ bp_range.cc error.o $(objs)

testArguments: arguments.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ arguments.cc error.o

testLogp: logp.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ logp.cc error.o $(objs)

testRand: rand.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ rand.cc error.o $(objs)

testsArray: sArray.cc sArray.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ sArray.cc error.o $(objs)

testcArray: cArray.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ cArray.cc error.o $(objs)


fileParserTest: fileParser.cc fileParser.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ fileParser.cc $(objs)

testHashDrivers: testHashMapList testHashTree testVHashSet testVHashMap testSHashSet testSHashMap testVSHashMap

testHashMapList: hash_map_list.cc hash_abstract.o hash_map_list.h hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ hash_map_list.cc error.o hash_abstract.o general.o

testHashTree: hash_tree.h hash_tree.cc hash_map_list.cc hash_abstract.o hash_map_list.h hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ hash_tree.cc error.o hash_abstract.o general.o

testVHashSet: vhash_set.h vhash_set.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ vhash_set.cc error.o hash_abstract.o general.o

testVHashMap: vhash_map.h vhash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ vhash_map.cc error.o hash_abstract.o general.o

testSHashSet: shash_set.h shash_set.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ shash_set.cc error.o hash_abstract.cc general.o 

testSHashMap: shash_map.h shash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ shash_map.cc error.o hash_abstract.cc general.o

testVSHashMap: vshash_map.h vshash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ vshash_map.cc error.o hash_abstract.cc general.o



#----------------------------------------------------------------------------
# Gather stuff together.
#----------------------------------------------------------------------------

all_srcs = error.cc bp_range.cc range.cc arguments.cc logp.cc rand.cc sArray.cc cArray.cc general.cc lineqsolve.cc debug.cc fileParser.cc hash_map_list.cc hash_abstract.cc  mul_mdmd_md.c mul_mdmd_md_l0g.c mul_mdmd_md_l0nf.c eig.c bin_search.cc # string_err.cc

all_hdrs = error.h debug.h fileParser.h hash_map_list.h hash_abstract.h eig.h

all_objs = 

all_exs = testError testDebug testBPRange testArguments testLogp testRand testsArray testcArray fileParserTest testHashDrivers testHashMapList testHashTree testVHashSet testVHashMap testSHashSet testSHashMap testVSHashMap



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
# when the file doesn't exist, do a 'touch depends.make' first.
depends: depend

# also include the dependencies
include depends.make


