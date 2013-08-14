
#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2013 Jeff Bilmes
# Licensed under the Open Software License version 3.0
#


# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://j.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([DBN.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.11 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

# Checks for programs.

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_AWK

# Checks for libraries.

m4_include([blas_args])
m4_include([blas_check])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

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
