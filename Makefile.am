EXTRA_DIST = RELEASE scripts
SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk tksrc

.PHONY: revdist datedist
revdist:
	make VERSION=`hg id | tr ' ' _` dist

datedist:
	make VERSION=`date +"%Y-%m-%d_%H%M"` dist
