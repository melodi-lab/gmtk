#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([adaptive_filter.h])
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
    AC_DEFINE([HAVE_HG_H],[1],[Distributed Mercurial revision id header availabl
e?])
  ],[
    AC_DEFINE([HGID],"UNKNOWN",[Unknown Mercurial revision id - broken distribut
ion?])
    AC_DEFINE([HGDATE],"UNKOWN",[Unkown Mercurial checkin date - broken distribution?])
  ])
])
AM_CONDITIONAL([HAVE_HG],[test  x$HAVE_HG = xyes])

# Checks for libraries.
AC_CHECK_LIB([m],[sqrt])

case "${host}" in
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRFTIME
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow sqrt strchr strerror strspn strtol strtoul log1p])

AC_ARG_WITH([logp], 
            [AS_HELP_STRING([--with-logp],
              [select logp implementation (log, table, smart, renorm) @<:@default: log@:>@])],
            [],
            [with_logp=log])
case "${with_logp}" in
log) 
  ;;
table)
  AC_DEFINE([_TABLE_], [1], [Define to use table logp implementation]) 
  ;;
smart)
  AC_MSG_ERROR([the smart logp table is not yet implemented]) 
  ;;
renorm)
  AC_MSG_ERROR([renormalization is not yet implemented]) 
  ;;
*)
  AC_MSG_ERROR([Must specify --with-logp= one of log, table, smart, renorm]) 
  ;;
esac

AC_ARG_WITH([LZERO],
            [AS_HELP_STRING([--with-LZERO],
              [the smallest log value that we represent @<:@default: -1.0E10@:>@])],
	    [],
            [with_LZERO=(-1.0E10)])

AC_DEFINE_UNQUOTED([LZERO], [(${with_LZERO})], [smallest log value that we represent])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
