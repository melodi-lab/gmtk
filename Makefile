#
# Top-level Makefile for GMTK
#
# $Header$
#


# other flags
EXLDFLAGS=
# GCC=/usr/nikola/pkgs/gcc/.3.2/bin/gcc
# GPP=/usr/nikola/pkgs/gcc/.3.2/bin/g++
CC=gcc
CXX=g++
OPTFLAGS = -g -O3 -Wno-deprecated -march=pentium3 -mfpmath=sse

# GMTK modules 
MODULES = \
	IEEEFloatingpoint \
	pfile \
	featureFileIO \
	miscSupport \
	tksrc

# files/dirs that should not be contained in distribution. 
EXCLUDE = \
	tests_and_development \
	EXCLUDE


MAKE_VARS = \
	OPTFLAGS="$(OPTFLAGS)" \
	CC="$(CC)" \
	CXX="$(CXX)" \
	EXLDFLAGS="$(EXLDFLAGS)"

all clean:
	for subdir in $(MODULES); do \
		(cd $$subdir; $(MAKE) $(MAKE_VARS) $@); \
	done

linux solaris ibm cygwin:
	for subdir in IEEEFloatingpoint; do \
		(cd $$subdir; $(MAKE) $(MAKE_VARS) $@); \
	done

depend:
	for subdir in $(MODULES); do \
		(cd $$subdir; touch depends.make; $(MAKE) $(MAKE_VARS) $@); \
	done


# this will make world just for linux for now.
World:
	$(MAKE) linux
	$(MAKE) clean
	$(MAKE) depend
	$(MAKE) all

TAR = /bin/tar

package:  EXCLUDE
	$(TAR) cvzXf EXCLUDE ../gmtk-`cat RELEASE`.tar.gz .

date:  EXCLUDE
	$(TAR) cvzXf EXCLUDE - . > ../gmtk-`date +%a_%b_%d_%Y_%k:%M | sed -e 's, ,,g'`.tar.gz

# always remake this target when called.
EXCLUDE: force
	(find $(EXCLUDE) -type d -print -prune ; \
	find . \( -name "*~" -o -name "*~[0-9]*" -o -name "core*" -o -name "*.o" -o -name "#*" -o -name ".#*" -o -name "*_bak" \) -print; \
	find $(MODULES) -type f -perm +0111 \! \( -name '*.cc' -o -name '*.h' \) ; \
	find . -name CVS -print; \
	find . -name RCS -print) | \
	sed 's,^\./,,' > $@



force:



#
#----------------------------------------------------------------------------
# CVS management.
#----------------------------------------------------------------------------

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
# Source Distribution Creation
#----------------------------------------------------------------------------

# not yet working.

srcdist:
#	/home/bilmes/bin/tar -cvf gmtkdist.tar `find $(src_dirs) -name bin -prune -o -name RCS -prune -o -name old -prune -o -type f \! \( -perm -u+x -o -perm -g+x -o -perm -o+x -o -name core -o -name \*.pfile -o -name \*.o -o -name \*.a -o -name \*.gmp -o -name \*~ -o -name gmtkEMtrain -o -name gmtkViterbi -o -name \*.pure -o -name \*.out -o -name \*,v -o -name \#\*\# \)`
	-@mv gmtkdist.tar.gz gmtkdist.previous.tar.gz 
	tar -cvf gmtkdist.tar Makefile gmtkmake `find $(src_dirs) -name bin -prune -o -name RCS -prune -o -name old -prune -o -name BUG\* -prune -o -name bug\* -prune -o -name exp\* -prune -o -name \*Tutorial\* -prune -o -type f \! \( -perm -u+x -o -perm -g+x -o -perm -o+x -o -name bug\* -o -name BUG\* -o -name core -o -name \*.pfile -o -name \*.o -o -name exp\* -o -name \*Tutorial\* -o -name \*.a -o -name \*.gmp -o -name \*~ -o -name gmtkEMtrain -o -name gmtkViterbi -o -name \*.pure -o -name \*.out -o -name \*,v -o -name \#\*\# -o -name \*.gz -o -name \*.uu -o -name \*.ps \) -print`
	gzip gmtkdist.tar


dist:
	-@mkdir $(DISTDIR)
	cp tksrc/gmtkEMtrain tksrc/gmtkViterbi tksrc/gmtkParmConvert tksrc/gmtkScore tksrc/gmtkSample $(DISTDIR)/.
	cp miscSupport/string_err $(DISTDIR)/.
	strip $(DISTDIR)/gmtk* $(DISTDIR)/string_err
	cp auroraTutorial/auroraTutorial.tar.gz $(DISTDIR)/.



# EoF