
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

EXTRA_DIST = RELEASE contrib version.sh version.info

dist_bin_SCRIPTS = \
scripts/fixTri.sh \
scripts/generate_random_graph.pl \
scripts/triangulateGA \
scripts/triangulateParallel \
scripts/triangulateTimings

SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk deepMLP tksrc tests

.PHONY: beta revdist datedist hgdatedist release beta-release

beta:
	make VERSION=`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr -cd A-Za-z0-9` dist

release: distcheck
	echo "PRETEND: hg commit -m \"This is the "`cat RELEASE`" release"
	echo "PRETEND: hg tag v"`cat RELEASE`
	echo "PRETEND: hg push"
	scp gmtk-`cat RELEASE`.tar.gz parrot:/g/melodi/html/web-pages/downloads/gmtk
	ssh parrot /g/melodi/html/web-pages/downloads/gmtk/mkhtml.sh

beta-release: distcheck
	echo "PRETEND: hg commit -m \"This is the "`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr -cd A-Za-z0-9`" release\""
	echo "PRETEND: hg tag v"`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr -cd A-Za-z0-9`
	echo "PRETEND: hg push"
	mv gmtk-`cat RELEASE`.tar.gz gmtk-`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr -cd A-Za-z0-9`.tar.gz
	scp gmtk-`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr -cd A-Za-z0-9`.tar.gz parrot:/g/melodi/html/web-pages/downloads/gmtk
	ssh parrot /g/melodi/html/web-pages/downloads/gmtk/mkhtml.sh

revdist:
	make VERSION=`cat RELEASE`+`hg id | tr ' ()' _` dist

datedist:
	make VERSION=`cat RELEASE`+`date +"%Y-%m-%d"` dist

hgdatedist:
	make VERSION=`cat RELEASE`+`date +"%Y-%m-%d"`.`hg id | tr ' ()' _` dist

