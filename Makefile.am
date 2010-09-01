# hgstamp will be created by the configure script
# It should contain the Mercurial revision id
# This line ensures it's included in the dist tarball
EXTRA_DIST = hgstamp

HG = @HG@
HGFLAGS = @HGFLAGS@
CMP = @CMP@
CMPFLAGS = @CMPFLAGS@

# this overwrites the revision id produced by $(HG)
# only if it differs from the current hgstamp.
# configure should take care that they never differ 
# in source trees unpacked from distributed tarballs
hgstamp: dummy
	test -f hgstamp || touch hgstamp
	$(HG) $(HGFLAGS) | $(CMP) $(CMPFLAGS) hgstamp - || \
	  $(HG) $(HGFLAGS) > hgstamp

.PHONY: dummy

SUBDIRS = IEEEFloatingpoint miscSupport featureFileIO mitk tksrc
