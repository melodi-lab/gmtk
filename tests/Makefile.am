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

LOCAL_GMTK_AT = $(wildcard gmtk_test_*.at)
EXTRA_DIST = testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in \
	     $(LOCAL_GMTK_AT)
TESTSUITE = testsuite

.PHONY: gmtk_tests.at
gmtk_tests.at:
	touch $@
	find $(srcdir) -type f -name gmtk_test_\*.at | \
	  $(AWK) '{print "m4_include(["$$1"])"}' > $@

check-local: atconfig atlocal $(TESTSUITE) gmtk_tests.at
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

installcheck-local: atconfig atlocal $(TESTSUITE) gmtk_tests.at
	$(SHELL) '$(TESTSUITE)' \
	AUTOTEST_PATH='$(bindir):$(abs_top_builddir)/tksrc' \
	$(TESTSUITEFLAGS)
#	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
#	$(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)' || \
	$(SHELL) '$(TESTSUITE)' --clean
	test ! -f test_and_dev.at || rm test_and_dev.at
	test ! -f gmtk_tests.at || rm gmtk_tests.at
	test ! -f $(TESTSUITE) || rm $(TESTSUITE)

distclean-local:
	test ! -f atconfig || rm atconfig

AUTOM4TE = $(SHELL) $(top_srcdir)/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/gmtk_tests.at test_and_dev.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@
