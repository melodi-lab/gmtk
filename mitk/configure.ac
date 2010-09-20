#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK], m4_esyscmd(echo -n `head -n 1 ../RELEASE`), [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CONFIG_SRCDIR([data-points.cc])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_CHECK_HEADERS([float.h math.h limits.h stdlib.h string.h sys/time.h values.h])

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
AC_CHECK_FUNCS([memset pow sqrt strchr strspn strtol])
# needed to make pfile.h happy
AC_FUNC_FSEEKO

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([discrete-mi])
AC_OUTPUT
