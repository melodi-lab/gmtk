#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# autoconf initialization
AC_PREREQ(2.59)
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}


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

# This identifies the CPU/OS to select the source files containing
# the correct platform-specific FPU code for the target. The source files
# are named ieeeFPsetup_@host_cpu@_@host_os@.{h,cc} in the Makefile.
AC_CANONICAL_HOST

# Turn on automake error checking and turn off whinging about missing
# some files GNU projects are supposed to have
AM_INIT_AUTOMAKE([1.11 silent-rules -Wall -Werror foreign])

# Check that we're in the correct directory by verifying the existence 
# of the file
AC_CONFIG_SRCDIR([ieeeFPsetup_ANY_ANY.cc])

# In case we need config.h to let the source code know about the 
# build environment
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
AC_PROG_RANLIB
AC_PROG_LN_S

AC_CHECK_FILE([../.hg/dirstate],[
  AC_CHECK_PROG([HAVE_HG],[hg],[yes])
  if test x$HAVE_HG = xyes; then
    HGREVID=`hg id ..`  
    HGDATE=`hg parents --template '{date|date}'`
    printf "#ifndef HGID\n#define HGID \"$HGREVID\"\n#define HGDATE \"$HGDATE\"\n#endif\n" > hgstamp.h
    AC_DEFINE([HAVE_HG_H],[1],[Live Mercurial revision id available?])
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

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([memset])

AM_CONDITIONAL([UNKNOWN_PLATFORM],[test ! -f "${srcdir}/ieeeFPsetup_${host_cpu}_${host_os}.cc"])

if test -f "${srcdir}/ieeeFPsetup_${host_cpu}_${host_os}.cc"; then
   echo Building IEEE floating point support for ${host_cpu}-${host_os}
else
   cp ${srcdir}/ieeeFPsetup_ANY_ANY.h \
      ieeeFPsetup_${host_cpu}_${host_os}.h
   cp ${srcdir}/ieeeFPsetup_ANY_ANY.cc \
      ieeeFPsetup_${host_cpu}_${host_os}.cc
   AC_WARN([
---------------------------------------------------------------
  Unable to catch internal floating point exceptions for ${host_cpu}-${host_os}
  GMTK will still build and run successfully, but with limited
  support for debugging floating point exceptions. Edit the 
  files IEEEFloatingpoint/ieeeFPsetup_${host_cpu}_${host_os}.{h,cc}
  to add reporting for floating point exceptions if you want
  to debug floating point exceptions on ${host_cpu}-${host_os}.
---------------------------------------------------------------
])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
