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

TESTSUITE = testsuite
LOCAL_GMTK_AT = \
gmtk_test_debug.at \
gmtk_test_newViterbi-1.at \
gmtk_test_newViterbi-2.at \
gmtk_test_newViterbi-3.at \
gmtk_test_newViterbi-4.at \
gmtk_test_ticket125.at \
gmtk_test_ticket127.at \
gmtk_test_ticket130.at \
gmtk_test_ticket136.at \
gmtk_test_ticket143.at \
gmtk_test_ticket149.at \
gmtk_test_ticket150.at \
gmtk_test_ticket162.at \
gmtk_test_ticket164.at \
gmtk_test_ticket166.at \
gmtk_test_ticket174.at \
gmtk_test_ticket191.at \
gmtk_test_ticket192.at \
gmtk_test_ticket213.at \
gmtk_test_ticket23.at \
gmtk_test_ticket310.at \
gmtk_test_ticket318.at \
gmtk_test_ticket320.at \
gmtk_test_ticket322.at \
gmtk_test_ticket326-1.at \
gmtk_test_ticket326-2.at \
gmtk_test_ticket326-3.at \
gmtk_test_ticket326-4.at \
gmtk_test_ticket348.at \
gmtk_test_ticket357.at \
gmtk_test_ticket42.at \
gmtk_test_ticket52.at

EXTRA_DIST = testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in \
	     gmtk_tests.at $(LOCAL_GMTK_AT)

.PHONY: gmtk_tests.at
gmtk_tests.at:
	touch $@
	find $(srcdir) -type f -name gmtk_test_\*.at | \
	  $(AWK) '{print "m4_include(["$$1"])"}' > $@

check-local: atconfig atlocal $(TESTSUITE) gmtk_tests.at
	$(SHELL) '$(srcdir)/$(TESTSUITE)' $(TESTSUITEFLAGS)

installcheck-local: atconfig atlocal $(TESTSUITE) gmtk_tests.at
	$(SHELL) '$(srcdir)/$(TESTSUITE)' \
	AUTOTEST_PATH='$(bindir):$(abs_top_builddir)/tksrc:$(abs_top_builddir)/miscSupport' \
	$(TESTSUITEFLAGS)

clean-local:
	test ! -f '$(TESTSUITE)' || \
	$(SHELL) '$(TESTSUITE)' --clean
	test ! -f test_and_dev.at || rm test_and_dev.at
	test ! -f gmtk_tests.at || rm gmtk_tests.at
	test ! -f $(TESTSUITE) || rm $(TESTSUITE)

distclean-local:
	test ! -f atconfig || rm atconfig
	test ! -f $(TESTSUITE).log || rm $(TESTSUITE).log 

AUTOM4TE = $(SHELL) $(top_srcdir)/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/gmtk_tests.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@
