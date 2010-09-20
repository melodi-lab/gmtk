#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([GMTK], m4_esyscmd(echo -n `head -n 1 ../RELEASE`), [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
AC_CONFIG_SRCDIR([pfile.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LEX
AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([float.h limits.h math.h stdlib.h string.h unistd.h inttypes.h stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_TYPE_OFF_T
AC_CHECK_SIZEOF([long])
AC_CHECK_DECLS([M_SQRT1_2,M_SQRT2],,,[#include <math.h>])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset sqrt strchr strcspn strerror strspn strstr strtol])
# does this go here?
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
