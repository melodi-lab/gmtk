
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2011 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

##### Autotest

# The `:;' works around a Bash 3.2 bug when the output is not writeable.
$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	:;{ \
	echo '# Signature of the current package.' && \
	echo 'm4_define([AT_PACKAGE_NAME],' && \
	echo '  [$(PACKAGE_NAME)])' && \
	echo 'm4_define([AT_PACKAGE_TARNAME],' && \
	echo '  [$(PACKAGE_TARNAME)])' && \
	echo 'm4_define([AT_PACKAGE_VERSION],' && \
	echo '  [$(PACKAGE_VERSION)])' && \
	echo 'm4_define([AT_PACKAGE_STRING],' && \
	echo '  [$(PACKAGE_STRING)])' && \
	echo 'm4_define([AT_PACKAGE_BUGREPORT],' && \
	echo '  [$(PACKAGE_BUGREPORT)])'; \
	echo 'm4_define([AT_PACKAGE_URL],' && \
	echo '  [$(PACKAGE_URL)])'; \
             } >'$(srcdir)/package.m4'

TESTSUITE = $(srcdir)/testsuite
TESTANDDEV= $(srcdir)/testanddev

include autoconf_tests

update-tests:
	(cd $(srcdir); ls gmtk_test_*.at) |                         \
	 $(AWK) 'BEGIN{print "LOCAL_GMTK_AT = \\"}                  \
	              {print $$0 " \\"}                             \
	           END{print "gmtk_tests.at"}' > autoconf_tests
	(cd $(srcdir); ls test_and_dev_scripts/test_and_dev_*.at) | \
	 $(AWK) 'BEGIN{print "LOCAL_TEST_AND_DEV = \\"}             \
	              {print $$0 " \\"}                             \
	           END{print "testanddev.at"}' >> autoconf_tests

EXTRA_DIST =          \
atlocal.in            \
testsuite.at          \
testanddev.at         \
$(srcdir)/package.m4  \
$(TESTSUITE)          \
$(TESTANDDEV)         \
$(LOCAL_GMTK_AT)      \
$(LOCAL_TEST_AND_DEV) \
autoconf_tests


$(srcdir)/gmtk_tests.at:
	touch $@
	find $(srcdir) -type f -name gmtk_test_\*.at | sort | \
	  $(AWK) '{print "m4_include(["$$1"])"}' > $@

$(srcdir)/testanddev.at:
	echo "AT_COPYRIGHT([Copyright (C) 2013 Jeff Bilmes])" > $@
	echo "AT_INIT([GMTK Test And Dev Suite])" >> $@
	echo "AT_COLOR_TESTS" >> $@
	find $(srcdir) -type f -name test_and_dev_\*.at | sort | \
	  $(AWK) '{print "m4_include(["$$1"])"}' >> $@

check-local: atconfig atlocal $(TESTSUITE) $(TESTANDDEV) gmtk_tests.at
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)


# The installcheck target is meant to test the installed build products
# rather than those sitting in the build tree. But the test suite uses
# some check_PROGRAMS that are not installed, so the test suite won't
# pass without the build tree in the AUTOTEST_PATH. We could create a
# separate test suite for the installed programs that doesn't require
# the extra test programs...

installcheck-local: atconfig atlocal $(TESTSUITE) $(TESTANDDEV) gmtk_tests.at
	$(SHELL) '$(TESTSUITE)' \
	AUTOTEST_PATH='$(bindir):$(abs_top_builddir)/tksrc:$(abs_top_builddir)/miscSupport:$(abs_top_builddir)/deepMLP' \
	$(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)'  || $(SHELL) '$(TESTSUITE)' --clean
	test ! -f test_and_dev.at || rm test_and_dev.at
	test ! -f gmtk_tests.at   || rm gmtk_tests.at
#	test ! -f $(TESTSUITE)    || rm $(TESTSUITE)
	test ! -f $(TESTANDDEV)   || $(SHELL) '$(TESTANDDEV)' --clean
	test ! -f testanddev.at   || rm testanddev.at
#	test ! -f $(TESTANDDEV)   || rm $(TESTANDDEV)

distclean-local:
	test ! -f atconfig || rm atconfig
	test ! -f $(TESTSUITE).log  || rm $(TESTSUITE).log 
	test ! -f $(TESTANDDEV).log || rm $(TESTANDDEV).log 

AUTOM4TE = $(SHELL) $(top_srcdir)/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/gmtk_tests.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@

$(TESTANDDEV): $(srcdir)/testanddev.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@
