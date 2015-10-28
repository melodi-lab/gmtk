
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2015 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

HELP2MANOPTS = -h -help -v -version  -N -I gmtk.h2m -n "a Graphical Models Toolkit (GMTK) "`cat $(srcdir)/../RELEASE`" program"
#--version-string=`cat $(srcdir)/../RELEASE`

SECTIONCMD = 's/\*\*\* \([A-Za-z0-9 ]*\) \*\*\*/.in -0.4i\n.B \1:\n.in 0.4i\n/g'

EXTRA_DIST = gmtk.h2m

dist_man1_MANS = \
  obs-print.1 \
  obs-info.1  \
  obs-skmeans.1 \
  obs-window.1 \
  obs-concat.1 \
  obs-diff.1 \
  obs-stats.1 \
  obs-cat.1 \
  gmtkTime.1 \
  gmtkJT.1 \
  gmtkEMtrain.1 \
  gmtkKernel.1 \
  gmtkViterbi.1 \
  gmtkTriangulate.1 \
  gmtkParmConvert.1 \
  gmtkTFmerge.1 \
  gmtkDTindex.1 \
  gmtkNGramIndex.1 \
  gmtkTie.1 \
  gmtkModelInfo.1 \
  gmtkOnline.1 \
  gmtkPrint.1 \
  gmtkMMItrain.1 \
  gmtkDMLPtrain.1 \
  gmtkViz.1

# Intentionally fails if you generate man pages in an 
# environment that doesn't build gmtkDMLPtrain or gmtkViz,
# because their man pages must exist for the distribution.

h2m_verbose   = $(h2m_verbose_@AM_V@)
h2m_verbose_  = $(h2m_verbose_@AM_DEFAULT_V@)
h2m_verbose_0 = @echo "  H2M     " $@;

obs-print.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-print$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-info.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-info$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-skmeans.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-skmeans$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-window.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-window$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-concat.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-concat$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-diff.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-diff$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-stats.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-stats$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-cat.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../featureFileIO/obs-cat$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@

gmtkTime.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkTime$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkJT.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkJT$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkEMtrain.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkEMtrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkKernel.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkKernel$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkViterbi.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkViterbi$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTriangulate.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkTriangulate$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkParmConvert.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkParmConvert$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTFmerge.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkTFmerge$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkDTindex.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkDTindex$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkNGramIndex.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkNGramIndex$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTie.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkTie$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkModelInfo.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkModelInfo$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkOnline.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkOnline$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkPrint.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkPrint$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkMMItrain.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkMMItrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkDMLPtrain.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkDMLPtrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkViz.1: gmtk.h2m
	$(h2m_verbose) $(HELP2MAN) $(HELP2MANOPTS) $(builddir)/../tksrc/gmtkViz$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
