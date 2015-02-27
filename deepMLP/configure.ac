
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2013 Jeff Bilmes
# Licensed under the Open Software License version 3.0
#


# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://j.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([DBN.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.12 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

AC_ARG_ENABLE([cppIfASCII],
              [AS_HELP_STRING([--disable-cppIfASCII],
                              [disable piping ASCII parameter files through cpp @<:@default: no@:>@])],
              [cppIfASCII=${enableval}], [cppIfASCII=yes])
if test x"$cppIfASCII" = x"yes"; then
  AC_DEFINE([PIPE_ASCII_FILES_THROUGH_CPP],[1],[Pipe ASCII files through cpp])
fi

# Checks for programs.

AC_PROG_CC
AC_PROG_CXX
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

AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_AWK

case "${host}" in
*cygwin*) AC_SUBST([XOPEN],["-D_XOPEN_SOURCE=900 -D__USE_XOPEN2K"]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" != x])


# Check for SSE headers
AC_CHECK_HEADERS([emmintrin.h], [have_sse=yes])

# Test if compiler accepts ISO C++ 2011.  This is required for
# the DMLP training code.

m4_include([m4/gmtk_cxx11])
if test x$gmtk_cxx11 = xyes -a x$have_sse = xyes -a x$XOPEN = x; then
  AC_SUBST([dmlp_library],[libDMLP.a])
  AC_SUBST([check_programs],['dmlp mnist label dmlpTest orderTest phipacTest'])
fi

# Checks for libraries.

m4_include([blas_args])
m4_include([blas_check])

# AX_PTHREAD comes from m4/m4_ax_pthread.m4 which is not distributed with GMTK
# See http://www.gnu.org/software/autoconf-archive/ax_pthread.html
AX_PTHREAD

# Checks for header files.

AC_CHECK_HEADERS([inttypes.h stdint.h])

# Checks for typedefs, structures, and compiler characteristics.

AC_TYPE_INT64_T

# Checks for library functions.

AC_CHECK_FUNCS([fdopen])
AC_FUNC_MMAP
AC_FUNC_FSEEKO
AC_SYS_LARGEFILE

AC_MSG_CHECKING([for cblas_daxpby])
  AC_LANG_PUSH([C++])
  mkl_save_CXXFLAGS=$CXXFLAGS
  CXXFLAGS=$BLAS_CFLAGS
  mkl_save_LIBS=$LIBS
  LIBS=$BLAS_LIBS

  if test "x$USE_MKL" = xyes; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "mkl.h"]],
         [[double x[3]; double y[3]; double alpha; double beta; cblas_daxpby(3, alpha, x, 1, beta, y, 1);]])],
      [AC_MSG_RESULT([yes])
       AC_DEFINE([HAVE_CBLAS_DAXPBY],1,[Define if CBLAS has cblas_daxpby])], 
      [AC_MSG_RESULT([no])])
  else
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[include <cblas.h>]],
         [[double x[3]; double y[3]; double alpha; double beta; cblas_daxpby(3, alpha, x, 1, beta, y, 1);]])],
      [AC_MSG_RESULT([yes])
       AC_DEFINE([HAVE_CBLAS_DAXPBY],1,[Define if CBLAS has cblas_daxpby])], 
      [AC_MSG_RESULT([no])])
  fi
  LIBS=$mkl_save_LIBS
  CXXFLAGS=$mkl_save_CXXFLAGS
  AC_LANG_POP([C++])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
