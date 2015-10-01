
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2015 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

HELP2MANOPTS = -h -help -v -version -N -I gmtk.h2m

SECTIONCMD = sed 's/\*\*\* \([A-Za-z0-9 ]*\) \*\*\*/.in -0.4i\n.B \1:\n.in 0.4i\n/g'

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
  $(dmlp_man1_mans) \
  $(gmtkViz_man1_mans)

EXTRA_man1_MANS = gmtkDMLPtrain.1 gmtkViz.1

obs-print.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-print$(EXEEXT)
obs-info.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-info$(EXEEXT)
obs-skmeans.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-skmeans$(EXEEXT)
obs-window.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-window$(EXEEXT)
obs-concat.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-concat$(EXEEXT)
obs-diff.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-diff$(EXEEXT)
obs-stats.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-stats$(EXEEXT)
obs-cat.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../featureFileIO/obs-cat$(EXEEXT)

gmtkTime.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkTime$(EXEEXT)
gmtkJT.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkJT$(EXEEXT)
gmtkEMtrain.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkEMtrain$(EXEEXT)
gmtkKernel.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkKernel$(EXEEXT)
gmtkViterbi.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkViterbi$(EXEEXT)
gmtkTriangulate.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkTriangulate$(EXEEXT)
gmtkParmConvert.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkParmConvert$(EXEEXT)
gmtkTFmerge.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkTFmerge$(EXEEXT)
gmtkDTindex.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkDTindex$(EXEEXT)
gmtkNGramIndex.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkNGramIndex$(EXEEXT)
gmtkTie.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkTie$(EXEEXT)
gmtkModelInfo.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkModelInfo$(EXEEXT)
gmtkOnline.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkOnline$(EXEEXT)
gmtkPrint.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkPrint$(EXEEXT)
gmtkMMItrain.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkMMItrain$(EXEEXT)
gmtkDMLPtrain.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkDMLPtrain$(EXEEXT)
gmtkViz.1: gmtk.h2m
	-$(HELP2MAN) --output=$@  $(HELP2MANOPTS) -n obs-cat $(builddir)/../tksrc/gmtkViz$(EXEEXT)
