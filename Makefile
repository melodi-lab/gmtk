#
# Top-level Makefile for GMTK
#
# $Header$
#

.EXPORT:
.EXPORT: nikola2 Linux

# other flags
EXLDFLAGS=
GCC=/usr/nikola/pkgs/gcc/.3.2/bin/gcc
GPP=/usr/nikola/pkgs/gcc/.3.2/bin/g++
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
	CC="$(GCC)" \
	CXX="$(GPP)" \
	EXLDFLAGS="$(EXLDFLAGS)"

depend all clean:
	for subdir in $(MODULES); do \
		(cd $$subdir; $(MAKE) $(MAKE_VARS) $@); \
	done

TAR = /bin/tar

package:  EXCLUDE
	$(TAR) cvzXf EXCLUDE ../gmtk-`cat RELEASE`.tar.gz .

EXCLUDE: force
	(find $(EXCLUDE) -type d -print -prune ; \
	find . \( -name "*~" -o -name "*~[0-9]*" -o -name "core*" -o -name "*.o" -o -name "#*" -o -name ".#*" \) -print; \
	find . -name CVS -print; \
	find . -name RCS -print) | \
	sed 's,^\./,,' > $@

force:


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