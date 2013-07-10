
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
#

EXTRA_DIST = RELEASE contrib 

dist_bin_SCRIPTS = \
scripts/fixTri.sh \
scripts/generate_random_graph.pl \
scripts/triangulateGA \
scripts/triangulateParallel \
scripts/triangulateTimings

SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk deepMLP tksrc tests

.PHONY: revdist datedist hgdatedist autoclean
revdist:
	make VERSION=`hg id | tr ' ()' _` dist

datedist:
	make VERSION=`date +"%Y-%m-%d_%H%M"` dist

hgdatedist:
	make VERSION=`date +"%Y-%m-%d"`.`hg id | tr ' ()' _` dist

