#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK],m4_esyscmd(echo -n `head -n 1 RELEASE`), [https://lungs.ee.washington.edu/trac/gmtk/newticket])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.

AC_CHECK_PROG([HAVE_HG],[hg],[yes])
AC_CHECK_FILE([.hg/dirstate],[
  # the builddir has the repository
  if test x"$HAVE_HG" = x"yes"; then
    # build host has Mercurial installed, so use it to keep
    # the revision id up to date
    AC_SUBST([HG],[hg])
    AC_SUBST([HGFLAGS],["parents --template '{node|short}'"])
    hg parents --template '{node|short}' > hgstamp 2> /dev/null || \
      echo Hg-ERROR > hgstamp
  else
    # how'd we get a repo w/o hg? use last id if we have it
    AC_CHECK_FILE([hgstamp],AC_SUBST([HG],[cat $(top_srcdir)/hgstamp]),
                            AC_SUBST([HG],[echo Hg-MISSING]))
  fi
],[
  # build dir does not have the repository, so use
  # distributed id if we have it
  AC_CHECK_FILE([hgstamp],AC_SUBST([HG],[cat $(top_srcdir)/hgstamp]),
                          AC_SUBST([HG],[echo ID-MISSING]))
])

AC_CHECK_PROG([HAVE_CMP],[cmp],[yes])
if test x"${HAVE_CMP}" = x"yes"; then
  AC_SUBST([CMP],[cmp])
  AC_SUBST([CMPFLAGS],[-s])
else
  AC_MSG_ERROR([cannot build GMTK without the cmp program])
fi

#AC_CHECK_PROG([HAVE_HG],[hg],[yes],[no])
#if test x"$HAVE_HG" = x"yes"; then
#  AC_CHECK_FILE([.hg/hgrc],[
#    # the build host has Mercurial installed and the 
#    # builddir has the repository, so use Mercurial
#    # to keep the revision id up to date
#    AC_SUBST([HG],[hg])
#    AC_SUBST([HGFLAGS],[parents --template '{node|short}'])
#    hg parents --template '{node|short}' > hgstamp 2> /dev/null || \
#      echo Hg-ERROR > hgstamp
#  ],[
#    # the build host has Mercurial but the builddir
#    # does not contain the repository, 
#    AC_SUBST([HG],[echo Hg-MISSING]
#    AC_SUBST([HGFLAGS],[])
#    echo Hg-MISSING > hgstamp]
#else 
#  AC_SUBST([HGFLAGS],[])
#  AC_CHECK_FILE([hgstamp],AC_SUBST([HG],[cat $(top_srcdir)/hgstamp],
#                          AC_SUBST([HG],[echo Hg-MISSING])
#fi

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([IEEEFloatingpoint
                   miscSupport
                   featureFileIO
                   mitk
                   tksrc])
AC_OUTPUT
