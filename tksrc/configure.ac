#                                               -*- Autoconf -*-

#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2010 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`), 
        [https://j.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([GMTK_DlinkMatrix.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.12 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

AC_ARG_ENABLE([register],
              [AS_HELP_STRING([--enable-register],
                              [use the deprecated register keyword @<:@default: no@:>@])],
              [register=${enableval}],[register=no])

if test x${register} = xno; then
  AC_DEFINE([REGISTER],[],[don't use the register keyword])
else
  AC_DEFINE([REGISTER],[register],[use the register keyword])
fi

AC_ARG_ENABLE([assertions],
              [AS_HELP_STRING([--disable-assertions],
                              [disable assertion checking @<:@default: no@:>@])],
              [assertions=${enableval}], [assertions=yes])
AC_ARG_ENABLE([gzip],
              [AS_HELP_STRING([--disable-gzip],
                              [disable piping .gz files through gzip @<:@default: no@:>@])],
              [gzip=${enableval}], [gzip=yes])
AC_ARG_ENABLE([bzip2],
              [AS_HELP_STRING([--disable-bzip2],
                              [disable piping .bz2 files through bzip2 @<:@default: no@:>@])],
              [bzip2=${enableval}], [bzip2=yes])
if test x"$gzip" = x"yes"; then
  AC_DEFINE([ENABLE_GZIP],[1],[Pipe .gz files through gzip])
fi
if test x"$bzip2" = x"yes"; then
  AC_DEFINE([ENABLE_BZIP2],[1],[Pipe .bz2 files through bzip2])
fi
AC_ARG_ENABLE([cppIfASCII],
	      [AS_HELP_STRING([--disable-cppIfASCII],
                              [disable piping ASCII parameter files through cpp @<:@default: no@:>@])],
              [cppIfASCII=${enableval}], [cppIfASCII=yes])
if test x"$cppIfASCII" = x"yes"; then
  AC_DEFINE([PIPE_ASCII_FILES_THROUGH_CPP],[1],[Pipe ASCII files through cpp])
fi

AC_ARG_ENABLE([tmpLocalPool],
              [AS_HELP_STRING([--disable-tmpLocalPool],
                              [disable use of temporary local clique value pool @<:@default: no@:>@])],
              [tmpLocalPool=${enableval}], [tmpLocalPool=yes])
if test x"$tmpLocalPool" = x"yes"; then
  AC_DEFINE([USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL],[1],[Hold clique values until after pruning])
fi

AC_ARG_ENABLE([continuous-cardinality-warning],
              [AS_HELP_STRING([--disable-continous-cardinality-warning],
                              [disable warning about continous variables with non-zero cardinality in trifiles written by older versions of GMTK @<:@default: no@:>@])],
              [continuous_cardinaliy_warning=${enableval}], [continuous_cardinality_warning=yes])
if test x"$continuous_cardinality_warning" = x"yes"; then
  AC_DEFINE([CONTINUOUS_CARDINALITY_WARNING],[1],[Warn about non-zero continuous RV cardinality in trifile])
fi

AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

AM_CONDITIONAL([NO_ASSERTIONS],[test x$assertions != xyes])
AM_CONDITIONAL([NO_MODEL_CHECKS],[test x$model_checks != xyes])

AC_ARG_WITH([LZERO],
            [AS_HELP_STRING([--with-LZERO],
              [the smallest log value that we represent (note that values > -1e16 may cause arithmetic problems (0+0 != 0) @<:@default: -1.0E17@:>@])],
	    [],
            [with_LZERO=-1.0E17])

badLZERO=`echo ${with_LZERO} | awk '{if ($1 > -1.0e16) {print "bad"}}'`
if test x"$badLZERO" = x"bad"; then
  AC_WARN([Values for LZERO greater than -1.0e16 may cause an arithmetic issue so that 0+0 != 0])
fi 
AC_DEFINE_UNQUOTED([LZERO], [(${with_LZERO})], [smallest log value that we represent])

m4_include([../deepMLP/blas_args])

# Checks for programs.

AC_PROG_CXX
AC_PROG_CC

#check for Clang, which lies about being GCC
AC_MSG_CHECKING([if compiling with clang])
AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([],[[
#ifndef __clang__
  not clang
#endif
]])],
[CLANG=yes],[CLANG=no])
AC_MSG_RESULT([$CLANG])
AM_CONDITIONAL([USE_GCC_FLAGS],[test x$GCC = xyes -a x$GXX = xyes -a x$CLANG = xno])
AM_CONDITIONAL([USE_CLANG_FLAGS],[test x$CLANG = xyes])

AM_PROG_LEX
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CHECK_FILE([$srcdir/../.hg/dirstate],[
  AC_CHECK_PROG([HAVE_HG],[hg],[yes])
  if test x$HAVE_HG = xyes; then
    HGREVID=`hg id $srcdir/..`  
    HGDATE=`hg parents --template '{date|date}'`
    printf "#ifndef HGID\n#define HGID \"$HGREVID\"\n#define HGDATE \"$HGDATE\"\n#endif\n" > hgstamp.h
    AC_DEFINE([HAVE_HG_H],[1],[Live Mercurial revision id header available?])
  fi
],[
  AC_CHECK_FILE([hgstamp.h],[
    AC_DEFINE([HAVE_HG_H],[1],[Distributed Mercurial revision id header available?])
  ],[
    AC_DEFINE([HGID],"UNKNOWN",[Unknown Mercurial revision id - broken distribution?])
    AC_DEFINE([HGDATE],"UNKOWN",[Unkown Mercurial checkin date - broken distribution?])
  ])
])
AM_CONDITIONAL([HAVE_HG],[test  x$HAVE_HG = xyes])

m4_include([m4/wxwin.m4])


AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])

if test x"$gmtkViz" = x"yes"; then
  WX_CONFIG_OPTIONS
  WX_CONFIG_CHECK([2.9.0], [HAVE_WXWIDGETS=yes], [AC_MSG_WARN([$wx_error_message

gmtkViz will not be built])])

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
fi
if test x"$HAVE_WXWIDGETS" = x"yes"; then
   AC_SUBST([GMTKVIZ],[gmtkViz])
   AC_SUBST([VIZCXXFLAGS],["$WX_CXXFLAGS $EXTRA_WX_FLAGS"])
   AC_SUBST([VIZLDFLAGS],[$WX_LIBS])
fi

case "${host}" in
*apple*)
  if test x"$HAVE_WXWIDGETS" = x"yes"; then
    AC_SUBST([OSXBUNDLE],[gmtkViz.app])

#   The EnableTop() bug (ticket 64) is fixed in wxWidgets 2.9.2
#   and gmtkViz souce checks if ENABLE_TOP_OFFSET is
#   needed using wxCHECK_VERSION(2,9,2)
    AC_DEFINE([GMTK_WX_OSX],[1],[Do we need to activate OSX fixes?])
  fi ;;
*-sun-solaris*) AC_CHECK_LIB([sunmath],[annuity]) ;;
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac

case "${host}" in
*apple*) ;; # don't use X11 on OS X
*) AC_SUBST([X_LIBS],[-lX11]) ;; 
esac

# Check for SSE headers
AC_CHECK_HEADERS([emmintrin.h], [have_sse=yes])

# Test if compiler accepts ISO C++ 2011.  This is required for
# the DMLP training code.

m4_include([../deepMLP/m4/gmtk_cxx11])
if test x$gmtk_cxx11 = xyes -a x$have_sse = xyes -a x$XOPEN = x; then
  AC_SUBST([dmlp_libraries],[../deepMLP/libDMLP.a])
  AC_SUBST([dmlp_programs],[gmtkDMLPtrain])
  AC_SUBST([dmlp_check_programs],["abstest dmlp"])
  AM_CONDITIONAL([BUILD_DMLP],[true])
else
  AM_CONDITIONAL([BUILD_DMLP],[false])
fi

# Checks for libraries.

# AX_PTHREAD comes from m4/m4_ax_pthread.m4 which is not distributed with GMTK
# See http://www.gnu.org/software/autoconf-archive/ax_pthread.html
AX_PTHREAD

     AC_ARG_WITH(
       [ps2eps],
       [AS_HELP_STRING([--with-ps2eps],
          [specify default ps2eps command ('no' or '--without-ps2eps' generates PostScript) @<:@default: ps2eps@:>@])],
       [AS_IF([test "x$with_ps2eps" != "xno"],
          [AC_DEFINE_UNQUOTED([PS2EPS_CMD], ["${with_ps2eps}"],[default ps2eps command])
           AC_DEFINE([EPS_EXT], ["eps"], [extension for print to file])
           AC_DEFINE([PRINT2FILE_ABBREV], ["EPS"], [abbreviation for print to file format])
           AC_DEFINE([PRINT2FILE_FORMAT], ["Encapsulated PostScript"], [format for print to file])
          ],
          [AC_DEFINE([PS2EPS_CMD], ["cat"], [pass-thru raw PostScript])
           AC_DEFINE([EPS_EXT], ["ps"], [extension for print to file])
           AC_DEFINE([PRINT2FILE_ABBREV], ["PS"], [abbreviation for print to file format])
           AC_DEFINE([PRINT2FILE_FORMAT], ["PostScript"], [format for print to file])
          ])
       ],
       [AC_DEFINE([PS2EPS_CMD], ["ps2eps"], [default ps2eps command])
        AC_DEFINE([EPS_EXT], ["eps"], [extension for print to file])
        AC_DEFINE([PRINT2FILE_ABBREV], ["EPS"], [abbreviation for print to file format])
        AC_DEFINE([PRINT2FILE_FORMAT], ["Encapsulated PostScript"], [format for print to file])
       ]
     )

AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])
AM_CONDITIONAL([MAKE_BUNDLE],[test  x$OSXBUNDLE = xgmtkViz.app])
AC_CHECK_LIB([m], [log])

m4_include([../m4/ax_lib_hdf5.m4])
AX_LIB_HDF5

AC_ARG_WITH([detBrush],
            [AS_HELP_STRING([--with-detBrush],
              [select gmtkViz brush for deterministic nodes (none, cross, grey) @<:@default=none@:>@])],
            [],
            [with_detBrush=none])
case "${with_detBrush}" in
none)
  AC_DEFINE([NO_DET_BRUSH],[1],[no special brush for deterministic variables])
  ;;
cross)
  AC_DEFINE([CROSS_DET_BRUSH],[1],[cross-hatch brush for deterministic variables])
  ;;
grey)
  AC_DEFINE([GREY_DET_BRUSH],[1],[solid grey brush for deterministic variables])
  ;;
*)
  AC_MSG_ERROR([Must specify --with-detBrush= one of none, cross, grey]) 
  ;;
esac

dnl The dlopen() function is in the C library for *BSD and in
dnl libdl on GLIBC-based systems
AC_SEARCH_LIBS([dlopen], [dl dld], [], [
  AC_MSG_ERROR([unable to find the dlopen() function])
])

m4_include([../deepMLP/blas_check])

# Checks for header files.
case "${host}" in
*apple*) ;; # don't use X11 on OS X
*) AC_PATH_XTRA ;;
esac

AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h dlfcn.h])

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
AC_CHECK_FUNCS([alarm floor memset pow regcomp sqrt strchr strerror strstr strtol log1p getline])

if test x$ac_cv_func_getline = x"yes"; then
  AC_LANG_PUSH([C++])

  AC_DEFUN([GMTK_GETLINE_PROG],[dnl
    AC_LANG_PROGRAM([dnl
#include <stdlib.h>
#include <stdio.h>
char *lineptr;
size_t n;
FILE *stream;
],[if (getline(&lineptr, &n, stream) > 0) exit(0)])])

  AC_MSG_CHECKING([if getline works without any feature test flags])
  AC_LINK_IFELSE([GMTK_GETLINE_PROG],[
    AC_MSG_RESULT([yes])
    AC_DEFINE([USE_GETLINE],[1],[use POSIX getline function])
  ],[
    AC_MSG_RESULT([no])
    AC_MSG_CHECKING([if getline requires _WITH_GETLINE])
    getline_check_save_CPPFLAGS=${CPPFLAGS}
    CPPFLAGS="-D_WITH_GETLINE ${getline_check_save_CPPFLAGS}"
    AC_LINK_IFELSE([GMTK_GETLINE_PROG],[
      AC_MSG_RESULT([yes])
      AC_DEFINE([_WITH_GETLINE],[1],[enable POSIX getline function])
      AC_DEFINE([USE_GETLINE],[1],[use POSIX getline function])
    ],[
      AC_MSG_RESULT([no])
      AC_MSG_CHECKING([if getline requires _POSIX_C_SOURCE >= 200809])
      CPPFLAGS="-D_POSIX_C_SOURCE=200809 ${getline_check_save_CPPFLAGS}"
      AC_LINK_IFELSE([GMTK_GETLINE_PROG],[
        AC_MSG_RESULT([yes])
        AC_DEFINE([_POSIX_C_SOURCE],[200809],[enable POSIX getline function])
        AC_DEFINE([USE_GETLINE],[1],[use POSIX getline function])
      ],[
        AC_MSG_RESULT([no])
        AC_MSG_WARN([unable to use getline function - will use fgets instead])
      ])
    ])
    CPPFLAGS=${getline_check_save_CPPFLAGS}
  ])
  AC_LANG_POP([C++])
fi

AC_FUNC_FSEEKO
AC_SYS_LARGEFILE

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
