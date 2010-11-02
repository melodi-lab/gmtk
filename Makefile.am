EXTRA_DIST = RELEASE contrib

dist_bin_SCRIPTS = \
scripts/generate_random_graph.pl \
scripts/triangulateGA \
scripts/triangulateParallel \
scripts/triangulateTimings

SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk tksrc

.PHONY: revdist datedist
revdist:
	make VERSION=`hg id | tr ' ' _` dist

datedist:
	make VERSION=`date +"%Y-%m-%d_%H%M"` dist
