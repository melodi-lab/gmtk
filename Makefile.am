
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

EXTRA_DIST = RELEASE contrib version.sh buildhelp.sh bootstrap.sh version.info gmtk.h2m

dist_bin_SCRIPTS = \
scripts/fixTri.sh \
scripts/generate_random_graph.pl \
scripts/triangulateGA \
scripts/triangulateParallel \
scripts/triangulateTimings

SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk deepMLP tksrc tests man

# machine to publish to (any that mount /g/melodi/html will probably do)
PUBLISHHOST=trumpet

# where to publish released tarballs
DOWNLOADS=$(PUBLISHHOST):/g/melodi/html/web-pages/downloads/gmtk

# script to update table of released tarballs
PUBLISH=/g/melodi/html/web-pages/downloads/gmtk/mkhtml.sh



.PHONY: beta revdist datedist hgdatedist release beta-release

# There should only be one release of a version number, so no -f on
# the hg tag. If it fails, you probably forgot to update the RELEASE file.
release: distcheck
	hg tag v`cat RELEASE` -m "This is the "`cat RELEASE`" release"
	hg push
	scp gmtk-`cat RELEASE`.tar.gz $(DOWNLOADS)
	ssh $(PUBLISHHOST) $(PUBLISH)

# Skip the distcheck because you've already done it and it's slow.
# You're on your honor not to skip it :)
quickrelease: dist
	hg tag v`cat RELEASE` -m "This is the "`cat RELEASE`" release"
	hg push
	scp gmtk-`cat RELEASE`.tar.gz $(DOWNLOADS)
	ssh $(PUBLISHHOST) $(PUBLISH)

# There can be lots of beta releases of the same version, even on the same
# date. So hg tag -f, and also include date and hg id in the tarball name.
beta-release: distcheck
	hg tag -f v`cat RELEASE`+`date +%Y%m%d` \
	  -m "This is the "`cat RELEASE`+`date +%Y%m%d`" release"
	hg push
	mv gmtk-`cat RELEASE`.tar.gz \
	  gmtk-`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr ' ()' _`.tar.gz
	scp gmtk-`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr ' ()' _`.tar.gz \
	  $(DOWNLOADS)
	ssh $(PUBLISHHOST) $(PUBLISH)

beta:
	make VERSION=`cat RELEASE`+`date +%Y%m%d`.`hg id -i | tr ' ()' _` dist

revdist:
	make VERSION=`cat RELEASE`+`hg id | tr ' ()' _` dist

datedist:
	make VERSION=`cat RELEASE`+`date +"%Y-%m-%d"` dist

hgdatedist:
	make VERSION=`cat RELEASE`+`date +"%Y-%m-%d"`.`hg id | tr ' ()' _` dist

