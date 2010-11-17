#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`), 
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([GMTK_DlinkMatrix.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LEX
AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CHECK_FILE([../.hg/dirstate],[
  AC_CHECK_PROG([HAVE_HG],[hg],[yes])
  if test x$HAVE_HG = xyes; then
    HGREVID=`hg id ..`  
    printf "#define HGID \"$HGREVID\"\n" > hgstamp.h
    AC_DEFINE([HAVE_HG_H],[1],[Live Mercurial revision id header available?])
  fi
],[
  AC_CHECK_FILE([hgstamp.h],[
    AC_DEFINE([HAVE_HG_H],[1],[Distributed Mercurial revision id header available?])
  ],[
    AC_DEFINE([HGID],"UNKNOWN",[Unknown Mercurial revision id - broken distribution?])
  ])
])
AM_CONDITIONAL([HAVE_HG],[test  x$HAVE_HG = xyes])

# Checks for libraries.

AC_CHECK_PROG([HAVE_WXCONFIG],[wx-config],[yes])
if test x"$HAVE_WXCONFIG" = x"yes"; then
  AC_LANG_PUSH([C++])
  wx_check_save_CXXFLAGS=${CXXFLAGS}
  CXXFLAGS=`wx-config --cxxflags`
  wx_check_save_LDFLAGS=${LDFLAGS}
  LDFLAGS=`wx-config --libs`
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include "wx/wx.h"],[wxString s("foo");s.Capitalize()])],
    AC_SUBST([HAVE_WXWIDGETS],[yes]),
    AC_MSG_WARN(
      [wx-config found but not wxWidgets >= 2.9. gmtkViz will not be built.])
  )
  LDFLAGS=${wx_check_save_LDFLAGS}
  CXXFLAGS=${wx_check_save_CXXFLAGS}
  AC_LANG_POP([C++])
fi
if test x"$HAVE_WXWIDGETS" = x"yes"; then
   AC_SUBST([GMTKVIZ],[gmtkViz])
   AC_SUBST([VIZCXXFLAGS],[`wx-config --cxxflags`])
   AC_SUBST([VIZLDFLAGS],[`wx-config --libs`])
fi

case "${host}" in
*apple*) 
  if test x"$HAVE_WXWIDGETS" = x"yes"; then
    AC_SUBST([OSXBUNDLE],[gmtkViz.app])
  fi ;;
*-sun-solaris*) AC_CHECK_LIB([sunmath],[annuity]) ;;
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])
AM_CONDITIONAL([MAKE_BUNDLE],[test  x$OSXBUNDLE = xgmtkViz.app])
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_PATH_XTRA
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_TYPE_OFF_T
AC_CHECK_SIZEOF([long])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([alarm floor memset pow regcomp sqrt strchr strerror strstr strtol])
AC_FUNC_FSEEKO

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
