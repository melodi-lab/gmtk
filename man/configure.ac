#                                               -*- Autoconf -*-

#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2015 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://j.ee.washington.edu/trac/gmtk/newticket])

AC_CONFIG_SRCDIR([gmtk.h2m])

AM_INIT_AUTOMAKE([1.12 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

# Checks for programs.

AC_PROG_CXX
AM_MISSING_PROG(HELP2MAN, help2man, $missing_dir)

#AC_PROG_INSTALL

# Test if compiler accepts ISO C++ 2011.  This is required for
# the DMLP training code.

m4_include([../deepMLP/m4/gmtk_cxx11])
if test x$gmtk_cxx11 = xyes -a x$have_sse = xyes -a x$XOPEN = x; then
  AC_SUBST([dmlp_man1_mans],[gmtkDMLPtrain.1])
fi

# Checks for libraries.

m4_include([../tksrc/m4/wxwin.m4])

WX_CONFIG_OPTIONS
WX_CONFIG_CHECK([2.9.0], [HAVE_WXWIDGETS=yes])

AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])

AC_DEFUN([GMTK_WCHAR_PROG],[dnl
AC_LANG_PROGRAM([dnl
#include <string>
using namespace std;
#include <wx/wx.h>
],[wxString s("foo")])])

if test x"$HAVE_WXWIDGETS" = x"yes" -a x$gmtkViz = xyes; then
  AC_MSG_CHECKING([whether wxWidgets works without -D_WCHAR_H_CPLUSPLUS_98_CONFORMANCE_])
  AC_LANG_PUSH([C++])
  wx_check_save_CXXFLAGS=${CXXFLAGS}
  wx_check_save_LIBS=${LIBS}
  CXXFLAGS="$wx_check_save_CXXFLAGS $WX_CXXFLAGS"
  LIBS="$wx_check_save_LIBS $WX_LIBS"
  AC_LINK_IFELSE([GMTK_WCHAR_PROG],[
    AC_MSG_RESULT([yes])
  ],[
    AC_MSG_RESULT([no])
    AC_MSG_CHECKING([whether wxWidgets works with -D_WCHAR_H_CPLUSPLUS_98_CONFORMANCE_])
    CXXFLAGS="$CXXFLAGS -D_WCHAR_H_CPLUSPLUS_98_CONFORMANCE_"
    AC_LINK_IFELSE([GMTK_WCHAR_PROG],[
      AC_MSG_RESULT([yes])
      EXTRA_WX_FLAGS=-D_WCHAR_H_CPLUSPLUS_98_CONFORMANCE_
    ],[ 
      AC_MSG_RESULT([no])
      AC_MSG_WARN([wxWidgets doesn't seem to work. gmtkViz will not be built.])
      HAVE_WXWIDGETS=no
    ])
  ])
  CXXFLAGS=$wx_check_save_CXXFLAGS
  LIBS=$wx_check_save_LIBS
  AC_LANG_POP([C++])
fi

if test x"$HAVE_WXWIDGETS" = x"yes"; then
   AC_SUBST([gmtkviz_man1_mans],[gmtkViz.1])
fi

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
