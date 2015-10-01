
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2015 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

HELP2MANOPTS = -h -help -v -version -N -I gmtk.h2m

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

obs-print.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-print $(builddir)/../featureFileIO/obs-print$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-info.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-info $(builddir)/../featureFileIO/obs-info$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-skmeans.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-skmeans $(builddir)/../featureFileIO/obs-skmeans$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-window.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-window $(builddir)/../featureFileIO/obs-window$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-concat.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-concat $(builddir)/../featureFileIO/obs-concat$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-diff.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-diff $(builddir)/../featureFileIO/obs-diff$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-stats.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-stats $(builddir)/../featureFileIO/obs-stats$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
obs-cat.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-cat$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@

gmtkTime.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkTime $(builddir)/../tksrc/gmtkTime$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkJT.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkJT $(builddir)/../tksrc/gmtkJT$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkEMtrain.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkEMtrain $(builddir)/../tksrc/gmtkEMtrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkKernel.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkKernel $(builddir)/../tksrc/gmtkKernel$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkViterbi.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkViterbi $(builddir)/../tksrc/gmtkViterbi$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTriangulate.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkTriangulate $(builddir)/../tksrc/gmtkTriangulate$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkParmConvert.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkParmConvert $(builddir)/../tksrc/gmtkParmConvert$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTFmerge.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkTFmerge $(builddir)/../tksrc/gmtkTFmerge$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkDTindex.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkDTindex $(builddir)/../tksrc/gmtkDTindex$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkNGramIndex.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkNGramIndex $(builddir)/../tksrc/gmtkNGramIndex$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkTie.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkTie $(builddir)/../tksrc/gmtkTie$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkModelInfo.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkModelInfo $(builddir)/../tksrc/gmtkModelInfo$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkOnline.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkOnline $(builddir)/../tksrc/gmtkOnline$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkPrint.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkPrint $(builddir)/../tksrc/gmtkPrint$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkMMItrain.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkMMItrain $(builddir)/../tksrc/gmtkMMItrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkDMLPtrain.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkDMLPtrain $(builddir)/../tksrc/gmtkDMLPtrain$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
gmtkViz.1: gmtk.h2m
	-$(HELP2MAN) $(HELP2MANOPTS) -n gmtkViz $(builddir)/../tksrc/gmtkViz$(EXEEXT) | \
	$(SED) $(SECTIONCMD) > $@
