#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([pfile.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([1.11 silent-rules -Wall -Werror foreign])

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
AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

AM_CONDITIONAL([NO_ASSERTIONS],[test x$assertions != xyes])
AM_CONDITIONAL([NO_MODEL_CHECKS],[test x$model_checks != xyes])
 
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LEX
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CHECK_FILE([../.hg/dirstate],[
  AC_CHECK_PROG([HAVE_HG],[hg],[yes])
  if test x$HAVE_HG = xyes; then
    HGREVID=`hg id ..`  
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


# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# this doesn't work - need to #include <H5Cpp.h>  and use namespace H5
#AC_LANG_PUSH([C++])
#AC_CHECK_LIB([hdf5_cpp], [H5File::isHdf5])  
#AC_LANG_POP([C++])

AC_MSG_CHECKING([for HDF5 >= 1.8])
AC_LANG_PUSH([C++])
hdf_check_save_LIBS=${LIBS}
LIBS="${LIBS} -lhdf5_cpp"
# H5Lregister should only be present in >= 1.8
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif],
[H5Lregister(NULL); H5File file("/dev/null", H5F_ACC_RDONLY )])],
        AC_SUBST([HAVE_HDF5],[yes]),[])
LIBS=${hdf_check_save_LIBS}
AC_LANG_POP([C++])
if test x"$HAVE_HDF5" = x"yes"; then 
  LIBS="-lhdf5_cpp $LIBS"
  AC_DEFINE([HAVE_LIBHDF5_CPP],[1],[Are the C++ bindings for HDF5 available?])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

case "${host}" in
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])

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
