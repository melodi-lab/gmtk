#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK], m4_esyscmd(echo -n `head -n 1 ../RELEASE`), [https://lungs.ee.washington.edu/trac/gmtk/newticket])
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

AC_CHECK_PROG([HAVE_WXWIDGETS],[wx-config],[yes],[no])
if test x"$HAVE_WXWIDGETS" = x"yes"; then
   AC_SUBST([GMTKVIZ],[gmtkViz])
   AC_SUBST([VIZCXXFLAGS],[`wx-config --cxxflags`])
   AC_SUBST([VIZLDFLAGS],[`wx-config --libs`])
fi

# Checks for libraries.
case "${host}" in
*-sun-solaris*) AC_CHECK_LIB([sunmath],[annuity]) ;;
esac
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_PATH_X
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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
