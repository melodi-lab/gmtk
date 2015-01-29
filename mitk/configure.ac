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

AC_CONFIG_SRCDIR([data-points.cc])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.12 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

AC_ARG_ENABLE([assertions],
              [AS_HELP_STRING([--disable-assertions],
                              [disable assertion checking @<:@default: no@:>@])],
              [assertions=${enableval}], [assertions=yes])

AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

AM_CONDITIONAL([NO_ASSERTIONS],[test x$assertions != xyes])
AM_CONDITIONAL([NO_MODEL_CHECKS],[test x$model_checks != xyes])

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

m4_ifdef([AM_PROG_AR],[AM_PROG_AR])

# Checks for libraries.
AC_CHECK_LIB([m], [log])

case "${host}" in
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])

# Checks for header files.
AC_CHECK_HEADERS([float.h math.h limits.h stdlib.h string.h sys/time.h values.h])

AC_LANG_PUSH([C++])
hdf_check_save_LIBS=${LIBS}
LIBS="${LIBS} -lhdf5_cpp -lhdf5"
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif],
[H5File file("/dev/null", H5F_ACC_RDONLY )])],
        AC_SUBST([HAVE_HDF5],[yes]),[])
LIBS=${hdf_check_save_LIBS}
AC_LANG_POP([C++])
if test x"$HAVE_HDF5" = x"yes"; then 
  LIBS="-lhdf5_cpp  -lhdf5 $LIBS"
  AC_DEFINE([HAVE_LIBHDF5_CPP],[1],[Are the C++ bindings for HDF5 available?])
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_TYPE_OFF_T
AC_CHECK_SIZEOF([long])
AC_CHECK_DECLS([DBL_MIN,MAXINT],[],[],[
#if HAVE_FLOAT_H
#include <float.h>
#endif
#if HAVE_LIMITS_H
#include <limits.h>
#endif
#if HAVE_VALUES_H
#include <values.h>
#endif
])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset pow sqrt finite strchr strspn strtol])
# needed to make pfile.h happy
AC_FUNC_FSEEKO

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([discrete-mi])
AC_OUTPUT
