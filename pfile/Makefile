
#change this to whatever
exec_prefix=./
docdir=./
topdir = ./
srcdir = ./


#these flags should ultimately be set by configure
#flags needed to set byte-swapping routines 
DFLAGS=-DINTV_WORDS_BIGENDIAN 
#limits definitions differ for solaris and linux
DFLAGS+= -DHAVE_IEEEFP_H
CPU=sparc

RANLIB = ranlib
INSTALL = /usr/bin/install -c
INSTALL_PROG = $(INSTALL)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_LIB = $(INSTALL_DATA)
MKINSTALLDIRS =	$(srcdir)/mkinstalldirs

#directories
INCLUDES= -I. -I../miscSupport -I../IEEEFloatingpoint
LIBDIR= -L. -L../miscSupport
PURIFYDIR=/g/ssli/commercial/rational/releases/DeveloperTools.DE17/sun4_solaris2/bin/

#new pfile library 
LIBS=-lm -lpfile  -lmiscSupport

#compilers and flags
CC=gcc
CXX=g++
PURIFY=purify
LD = $(CXX)
# CFLAGS = -g -O 
# CXXFLAGS = -g -O 
OPTFLAGS=-O3


CFLAGS = -g $(OPTFLAGS)
CXXFLAGS = -g $(OPTFLAGS)
LDFLAGS = $(LIBDIR)
#PUR=$(PURIFYDIR)/$(PURIFY)
PUR=

compile.c =  $(PUR) $(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES)
compile.cxx =  $(PUR) $(CXX) $(CXXFLAGS) $(DFLAGS) $(INCLUDES)

pure = 
link = $(PUR) $(LD) $(LDFLAGS) $(CXXFLAGS)

.c.o:
	$(compile.c) -c $<

.C.o:
	$(compile.cxx) -c $<

.cc.o:
	$(compile.cxx) -c $<


#----------------------------------------------------------------------------
# Default target.
#----------------------------------------------------------------------------

all: libpfile.a pfile_gaussian pfile_info pfile_merge pfile_diff pfile_print \
	pfile_select pfile_stats pfile_build  pfile_create \
	pfile_concat pfile_norm pfile_normutts pfile_ftrcombo  \
	pfile_klt pfile_mi pfile_skmeans pfile_addsil pfile_initmg


#----------------------------------------------------------------------------
# Build pfile_utils programs.
#----------------------------------------------------------------------------

libpfile.a: pfile.o spi.o Range.o error.o 
	ar rv libpfile.a pfile.o spi.o Range.o 
	-ranlib libpfile.a

pfile_skmeans : pfile_skmeans.o Range.o libpfile.a
	$(link) -o pfile_skmeans pfile_skmeans.o Range.o rand.o $(LIBS)

pfile_gaussian : pfile_gaussian.o bin_search.o libpfile.a
	$(link) -o pfile_gaussian pfile_gaussian.o bin_search.o $(LIBS)

#pfile_initmg : pfile_initmg.o rand.o libpfile.a
#	$(link) -o pfile_initmg pfile_initmg.o rand.o $(LIBS)

pfile_info : pfile_info.o libpfile.a
	$(link) -o pfile_info pfile_info.o $(LIBS)


pfile_initmg : pfile_initmg.o libpfile.a
	$(link) -o pfile_initmg pfile_initmg.o $(LIBS)




pfile_merge : pfile_merge.o libpfile.a
	$(link) -o pfile_merge pfile_merge.o $(LIBS)

pfile_ftrcombo : pfile_ftrcombo.o libpfile.a
	$(link) -o pfile_ftrcombo pfile_ftrcombo.o $(LIBS)

pfile_diff : pfile_diff.o libpfile.a
	$(link) -o pfile_diff pfile_diff.o $(LIBS)

pfile_print : pfile_print.o libpfile.a
	$(link) -o pfile_print pfile_print.o $(LIBS)

pfile_process : pfile_process.o libpfile.a
	$(link) -o pfile_process pfile_process.o $(LIBS)

pfile_select : pfile_select.o libpfile.a
	$(link) -o pfile_select pfile_select.o $(LIBS)

pfile_addsil : pfile_addsil.o rand.o
	$(link) -o pfile_addsil pfile_addsil.o rand.o $(LIBS)

pfile_norm : pfile_norm.o
	$(link) -o pfile_norm pfile_norm.o $(LIBS)

pfile_normutts : pfile_normutts.o 
	$(link) -o pfile_normutts pfile_normutts.o $(LIBS)

pfile_stats : pfile_stats.o
	$(link) -o pfile_stats pfile_stats.o  $(LIBS)

pfile_build : pfile_build.o
	$(link) -o pfile_build pfile_build.o $(LIBS)

pfile_create : pfile_create.o
	$(link) -o pfile_create pfile_create.o $(LIBS)

pfile_concat : pfile_concat.o  
	$(link) -o pfile_concat pfile_concat.o  $(LIBS)

pfile_addframes : pfile_addframes.o
	$(link) -o pfile_addframes pfile_addframes.o  $(LIBS)

pfile2rapbin : pfile2rapbin.o 
	$(link) -o pfile2rapbin pfile2rapbin.o $(LIBS)

pfile_eval : pfile_eval.o
	$(link) -o pfile_eval pfile_eval.o  $(LIBS)

pfile_combine : pfile_combine.o 
	$(link) -o pfile_combine pfile_combine.o $(LIBS)
pfile_corr : pfile_corr.o 
	$(link) -o pfile_corr pfile_corr.o $(LIBS)
pfile_subsample : pfile_subsample.o 
	$(link) -o pfile_subsample pfile_subsample.o $(LIBS)
pfile_mi : pfile_mi.o  Range.o pfile.o  MixBiNormal.o MixBiNormal_chelp.o rand.o spi.o
	$(link) -o pfile_mi pfile_mi.o Range.o pfile.o  MixBiNormal.o MixBiNormal_chelp.o rand.o spi.o  $(LIBS)

MixBiNormal: $(srcdir)/MixBiNormal.cc MixBiNormal_chelp.o rand.o
	$(compile.cxx) -DMAIN -o MixBiNormal $(srcdir)/MixBiNormal.cc pfile.o  MixBiNormal_chelp.o rand.o $(LIBS)


pfile_patchlabels : pfile_patchlabels.o  pfile.o 
	$(link) -o pfile_patchlabels pfile_patchlabels.o pfile.o  $(LIBS)

pfile_klt : pfile_klt.o Range.o pfile.o  mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o eig.o
	$(link) -o pfile_klt pfile_klt.o Range.o pfile.o  mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o eig.o $(LIBS)


#----------------------------------------------------------------------------
# Gather stuff together.
#----------------------------------------------------------------------------

all_srcs = $(srcdir)/bin_search.cc $(srcdir)/parse_subset.cc $(srcdir)/pfile_gaussian.cc $(srcdir)/pfile_info.cc $(srcdir)/pfile_merge.cc $(srcdir)/pfile_ftrcombo.cc $(srcdir)/pfile_diff.cc $(srcdir)/pfile_print.cc $(srcdir)/pfile_process.cc  $(srcdir)/pfile_select.cc $(srcdir)/pfile_norm.cc $(srcdir)/pfile_normutts.cc $(srcdir)/pfile_stats.cc $(srcdir)/pfile_mi.cc  $(srcdir)/pfile_build.cc $(srcdir)/pfile_create.cc  $(srcdir)/pfile_concat.cc $(srcdir)/pfile_patchlabels.cc  $(srcdir)/pfile_klt.cc $(srcdir)/pfile_skmeans.cc $(srcdir)/pfile_addsil.cc $(srcdir)/mul_mdmd_md.c $(srcdir)/mul_mdmd_md_l0g.c  $(srcdir)/mul_mdmd_md_l0nf.c $(srcdir)/eig.c $(srcdir)/Range.C $(srcdir)/rand.cc $(srcdir)/MixBiNormal.cc $(srcdir)/MixBiNormal_chelp.c $(srcdir)/pfile_initmg.cc


all_hdrs = $(srcdir)/eig.h  $(srcdir)/parse_subset.h $(srcdir)/Range.H $(srcdir)/rand.h $(srcdir)/MixBiNormal.h

all_objs = bin_search.o parse_subset.o pfile_gaussian.o pfile_info.o pfile_merge.o pfile_ftrcombo.o pfile_diff.o pfile_print.o pfile_process.o pfile_select.o pfile_norm.o pfile_normutts.o pfile_stats.o pfile_mi.o pfile_build.o pfile_create.o  pfile_concat.o pfile_addsil.o pfile_patchlabels.o pfile_klt.o pfile_skmeans.o pfile.o   mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o eig.o Range.o MixBiNormal.o MixBiNormal_chelp.o rand.o pfile_initmg.o pfile.o spi.o


# all_exs = pfile_gaussian pfile_info pfile_merge pfile_ftrcombo pfile_diff pfile_print pfile_process pfile_select pfile_norm pfile_normutts pfile_stats pfile_mi pfile_build pfile_klt pfile_create ilab_create pfile_skmeans pfile_concat pfile_addsil MixBiNormal pfile_initmg

all_exs = pfile_gaussian pfile_info pfile_merge pfile_ftrcombo pfile_diff pfile_print pfile_process pfile_select pfile_norm pfile_normutts pfile_stats pfile_mi pfile_build pfile_klt pfile_create  pfile_skmeans pfile_concat pfile_addsil pfile_initmg

all_docs = $(srcdir)/pfile_utils.1 README TODO

#----------------------------------------------------------------------------
# Installation
#----------------------------------------------------------------------------

install-docs:
	$(INSTALL_DATA) $(srcdir)/pfile_utils.1 $(docdir)/pfile_utils.1
	(cd $(docdir); \
	 ln -f pfile_utils.1 pfile_gaussian.1; \
	 ln -f pfile_utils.1 pfile_info.1; \
	 ln -f pfile_utils.1 pfile_merge.1; \
	 ln -f pfile_utils.1 pfile_ftrcombo.1; \
	 ln -f pfile_utils.1 pfile_concat.1; \
	 ln -f pfile_utils.1 pfile_diff.1; \
	 ln -f pfile_utils.1 pfile_print.1; \
	 ln -f pfile_utils.1 pfile_select.1; \
	 ln -f pfile_utils.1 pfile_norm.1; \
	 ln -f pfile_utils.1 pfile_normutts.1; \
	 ln -f pfile_utils.1 pfile_stats.1; \
	 ln -f pfile_utils.1 pfile_mi.1; \
	 ln -f pfile_utils.1 pfile_build.1; \
	 ln -f pfile_utils.1 pfile_patchlabels.1; \
	 ln -f pfile_utils.1 pfile_klt.1;)


uninstall-docs:
	rm -f $(docdir)/pfile_gaussian.1 $(docdir)/pfile_info.1 \
		$(docdir)/pfile_merge.1 $(docdir)/pfile_ftrcombo.1 \
		$(docdir)/pfile_concat.1 $(docdir)/pfile_diff.1 \
		$(docdir)/pfile_print.1 $(docdir)/pfile_select.1 \
		$(docdir)/pfile_norm.1 $(docdir)/pfile_normutts.1 \
		$(docdir)/pfile_stats.1 $(docdir)/pfile_mi.1 \
		$(docdir)/pfile_build.1 $(docdir)/pfile_patchlabels.1 \
		$(docdir)/pfile_klt.1

install: $(all_exs) install-docs
	$(INSTALL_PROG) pfile_gaussian $(exec_prefix)/bin/pfile_gaussian
	$(INSTALL_PROG) pfile_info $(exec_prefix)/bin/pfile_info
	$(INSTALL_PROG) pfile_merge $(exec_prefix)/bin/pfile_merge
	$(INSTALL_PROG) pfile_ftrcombo $(exec_prefix)/bin/pfile_ftrcombo
	$(INSTALL_PROG) pfile_diff $(exec_prefix)/bin/pfile_diff
	$(INSTALL_PROG) pfile_print $(exec_prefix)/bin/pfile_print
	$(INSTALL_PROG) pfile_select $(exec_prefix)/bin/pfile_select
	$(INSTALL_PROG) pfile_norm $(exec_prefix)/bin/pfile_norm
	$(INSTALL_PROG) pfile_normutts $(exec_prefix)/bin/pfile_normutts
	$(INSTALL_PROG) pfile_stats $(exec_prefix)/bin/pfile_stats
	$(INSTALL_PROG) pfile_mi $(exec_prefix)/bin/pfile_mi
	$(INSTALL_PROG) pfile_skmeans $(exec_prefix)/bin/pfile_skmeans
	$(INSTALL_PROG) pfile_addsil $(exec_prefix)/bin/pfile_addsil
	$(INSTALL_PROG) pfile_build $(exec_prefix)/bin/pfile_build
	$(INSTALL_PROG) pfile_create $(exec_prefix)/bin/pfile_create
	$(INSTALL_PROG) pfile_concat $(exec_prefix)/bin/pfile_concat
#	$(INSTALL_PROG) pfile_patchlabels $(exec_prefix)/bin/pfile_patchlabels
	$(INSTALL_PROG) pfile_klt $(exec_prefix)/bin/pfile_klt
	$(INSTALL_PROG) pfile_initmg $(exec_prefix)/bin/pfile_initmg

uninstall: uninstall-docs
	rm -f $(exec_prefix)/bin/pfile_gaussian
	rm -f $(exec_prefix)/bin/pfile_info
	rm -f $(exec_prefix)/bin/pfile_merge
	rm -f $(exec_prefix)/bin/pfile_ftrcombo
	rm -f $(exec_prefix)/bin/pfile_diff
	rm -f $(exec_prefix)/bin/pfile_print
	rm -f $(exec_prefix)/bin/pfile_select
	rm -f $(exec_prefix)/bin/pfile_norm
	rm -f $(exec_prefix)/bin/pfile_normutts
	rm -f $(exec_prefix)/bin/pfile_stats
	rm -f $(exec_prefix)/bin/pfile_mi
	rm -f $(exec_prefix)/bin/pfile_skmeans
	rm -f $(exec_prefix)/bin/pfile_addsil
	rm -f $(exec_prefix)/bin/pfile_build
	rm -f $(exec_prefix)/bin/pfile_create
	rm -f $(exec_prefix)/bin/pfile_concat
	rm -f $(exec_prefix)/bin/pfile_patchlabels
	rm -f $(exec_prefix)/bin/pfile_klt


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

rcs_files = $(all_srcs) $(all_hdrs) $(doc_files) $(config_srcs)

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

lib_name=libpfile.a
target_junk = $(all_objs) $(all_exs) $(lib_name)

clean:
	rm -rf $(target_junk)

distclean: clean
	rm -f config.log config.status config.cache Makefile

maintainer_junk = $(srcdir)/TAGS

maintainer-clean: distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."
	rm -rf $(maintainer_junk)

#################################################################
## Package distribution.
#################################################################

VERSION = 0.1

# The name of the distribution
distname = pfile_utils-$(VERSION)

# The file name of the compressed tar file
distfile = $(distname).tar.gz

# Where to build the disttree
distroot = .

# The directory in which the distribution is built
distdir = $(distroot)/$(distname)

# Build a distribution
dist: $(distfile)

# Build the file that contains the distribution
$(distfile): disttree
	rm -f $(distfile); \
	(cd $(distroot); tar cvf - $(distname) ) | gzip >$(distfile)

dist_files =	$(all_srcs) $(all_hdrs) $(all_docs) \
		$(config_srcs) $(config_utils) $(config_outfiles)

# Build a directory tree that can be tared and gzipped to form the distribution
disttree: $(dist_files)
	rm -fr $(distdir); \
	mkdir $(distdir); \
	chmod a=rwx $(distdir); \
	for i in $(dist_files) ; do \
		file=`echo $$i | sed -e "s^$(srcdir)/^^"`; \
		cp -p $(srcdir)/$$file $(distdir)/$$file || exit 1; \
	done;

# DO NOT DELETE THIS LINE -- make depend depends on it.
