#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(GMTK, 0.0, https://lungs.ee.washington.edu/trac/gmtk/newticket)
: ${CXXFLAGS=""}
AC_CONFIG_SRCDIR([mixNormalCollection.cc])
AC_CONFIG_HEADER([ac_config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h values.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_SIZEOF([long])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([strchr strtol])

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
