#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(GMTK, 0.0, https://lungs.ee.washington.edu/trac/gmtk/newticket)
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
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
#AC_TYPE_INT32_T
#AC_TYPE_UINT32_T
#AC_TYPE_INT64_T
#AC_TYPE_LONG_LONG_INT
AC_CHECK_SIZEOF([long])

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
