
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


# We need dgemm, which can come from PHiPAC, MKL, or BLAS
# PHiPAC is "built-in" to GMTK, so it is preferred. A well-tuned
# MKL or BLAS dgemm might be faster, so --disable-phipac if
# you don't want the PHiPAC dgemm

AC_ARG_ENABLE([phipac],
   [AS_HELP_STRING([--disable-phipac], [use BLAS or MKL matrix multiply instead of PHiPAC])],
   [AS_IF([test x$enable_phipac == xyes], 
      [AC_DEFINE([USE_PHIPAC],1,[prefer PHiPAC dGEMM over MKL/BLAS])])],
   [AC_DEFINE([USE_PHIPAC],1,[prefer PHiPAC dGEMM over MKL/BLAS])
    PHiPAC_LIB=libPHiPAC.a
    AC_SUBST(PHiPAC_LIB)])

# We need a BLAS. ATLAS or Goto are nice. MKL also provides BLAS.
# By default, we'll try to use a BLAS located either by pkg-config
# or in your default build search paths. You can specify
#   --with-blas-include=<dir> 
#   --with-blas-libdir=<dir>
#   --with-blas-libs=<libraries>
# if you want to use a particular BLAS. If you don't specify any
# of the above, we will try pkg-config cblas to locate the BLAS.
# If the --with-blas-* weren't specified and pkg-config fails to 
# find the BLAS, we just guess -I/usr/include -L/usr/lib -lcblas
#
# If you prefer to use MKL instead of a BLAS, see --with-mkl

AC_ARG_WITH([blas-include],
   [AS_HELP_STRING([--with-blas-include=<dir>], [use non-MKL BLAS includes found under <dir> @<:@/usr/include@:>@])],
   [AS_IF([test "x$with_blas_include" != "xno"], [BLAS_INCLUDE=-I$with_blas_include], [BLAS_INCLUDE=""])
    BLAS_REQUEST=yes],
   [BLAS_INCLUDE=""])
# assume -I/usr/include is already included in the include seach path  

AC_ARG_WITH([blas-libdir],
   [AS_HELP_STRING([--with-blas-libdir=<dir>], [use non-MKL BLAS libraries found under <dir> @<:@/usr/lib@:>@])],
   [AS_IF([test "x$with_blas_libdir" != "xno"], [BLAS_LIBDIR=-L$with_blas_libdir], [BLAS_LIBDIR=""])
    BLAS_REQUEST=yes],
   [BLAS_LIBDIR=""])
# assume -L/usr/lib is already included in the library seach path  

AC_ARG_WITH([blas-libs],
   [AS_HELP_STRING([--with-blas-libs=<libraries>], [-l flags for non-MKL BLAS library @<:@-lcblas@:>@])],
   [AS_IF([test "x$with_blas_libs" != "xno"], [BLAS_LIBS=$with_blas_libs], [BLAS_LIBS=""])
    BLAS_REQUEST=yes],
   [BLAS_LIBS=-lcblas])
# assume libcblas.* is what we want


#AC_ARG_WITH([blas],
#   [AS_HELP_STRING([--with-blas=<dir>], [use (non-MKL) BLAS installation in <dir>])],
#   [AS_IF([test "x$with_blas" != "xno"],[use_blas=yes; BLASROOT=$with_blas],[use_blas=no])],
#   [use_blas=no])

# If you want to use MKL instead of BLAS, use --with-mkl=<dir> to
# specify where the MKL installation lives. 

# NOTE that it is the absence of --with-mk (or presence of --without-mkl)
# that turns on usage of BLAS. 

AC_ARG_WITH([mkl],
   [AS_HELP_STRING([--with-mkl=<dir>], [use MKL installation in <dir>])],
   [AS_IF([test "x$with_mkl" != "xno"],[use_mkl=yes; MKLROOT=$with_mkl; use_blas=no],[use_mkl=no; use_blas=yes])],
   [use_mkl=no; use_blas=yes])

AC_ARG_WITH([mkl-thread],
   [AS_HELP_STRING([--with-mkl-thread=<lib>], 
     [use MKL threading library <lib> (sequential, intel_thread, pgi_thread, gnu_thread) @<:@sequential@:>@])],
   [case $with_mkl_thread in
      intel_thread | pgi_thread | sequential | gnu_thread) ;;
      no) with_mkl_thread=sequential ;;
      *) AC_MSG_ERROR([Unknown MKL threading layer: $with_mkl_threading]) ;;
    esac],
    [with_mkl_thread=sequential])

AC_ARG_ENABLE([mkl-ilp64],
  [AS_HELP_STRING([--enable-mkl-ilp64], [use ilp64 instead of lp64 interface for MKL])],
  [AS_IF([test "x$enable_mkl_ilp64" != "xno"], [mkl_interface=ilp64], [mkl_interface=lp64])],
  [mkl_interface=lp64])

AC_ARG_ENABLE([mkl-ia32],
  [AS_HELP_STRING([--enable-mkl-ia32], [use ia32 instead of x84_64 build of MKL])],
  [AS_IF([test "x$enable_mkl_ia32" != "xno"], [mkl_arch=ia32], [mkl_arch=intel64])],
  [mkl_arch=intel64])

AC_ARG_WITH([boost-include],
  [AS_HELP_STRING([--with-boost-include=<dir>], [use Boost headers in <dir> @<:@/usr/include@:>@])],
  [AS_IF([test x$with_boot_include != xno],[BOOST_CPPFLAGS=-I$with_boost_include])])

AC_ARG_WITH([boost-libdir],
  [AS_HELP_STRING([--with-boost-libdir=<dir>], [use Boost libraries in <dir> @<:@/usr/lib@:>@])],
  [AS_IF([test x$with_boot_libdir != xno],[BOOST_LDFLAGS=-L$with_boost_libdir])])

AC_SUBST([BOOST_CPPFLAGS])
AC_SUBST([BOOST_LDFLAGS])

# Checks for programs.

AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_AWK

# Checks for libraries.

if test x"$use_mkl" = xyes; then     # --with-mkl=<dir>
  if test x"$use_blas" = xyes; then
    AC_MSG_ERROR([Cannot use both MKL and BLAS])  # configure script logic should make this impossible
  fi
  BLAS_LIBS="-Wl,--start-group -lmkl_intel_${mkl_interface} -lmkl_${with_mkl_thread} -lmkl_core -Wl,--end-group"
  BLAS_CPPFLAGS=-I${MKLROOT}/include
  BLAS_LDFLAGS=-L${MKLROOT}/lib/${mkl_arch}

  AC_MSG_CHECKING([for MKL under $MKLROOT])
  AC_LANG_PUSH([C++])
  mkl_save_CXXFLAGS=$CXXFLAGS
  CXXFLAGS=$BLAS_CPPFLAGS
  mkl_save_LIBS=$LIBS
  LIBS="$BLAS_LDFLAGS $BLAS_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "mkl.h"]],
                                  [[double x[3]; double y[3]; double alpha; cblas_daxpy(3, alpha, x, 1, y, 1);]])],
    [AC_MSG_RESULT([yes])
     AC_DEFINE([HAVE_MKL],1,[Define if you have intel MKL])], 
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([MKL not found under $MKLROOT])])
  LIBS=$mkl_save_LIBS
  CXXFLAGS=$mkl_save_CXXFLAGS
  AC_LANG_POP([C++])
else
  # if any --with-blas-* check @ user specified location
  if test x$BLAS_REQUEST == xyes; then
     AC_MSG_CHECKING([for BLAS])
     AC_LANG_PUSH([C])
     blas_save_CFLAGS=${CFLAGS}
     CFLAGS=$BLAS_INCLUDE
     blas_save_LIBS=${LIBS}
     LIBS="$BLAS_LIBDIR $BLAS_LIBS"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                       [[double x[3]; double y[3]; double alpha; cblas_daxpy(3, alpha, x, 1, y, 1);]])],
       [AC_MSG_RESULT([yes])
        AC_DEFINE([HAVE_BLAS],1,[Define if you have a (non-MKL) BLAS installation])
        BLAS_CPPFLAGS=$BLAS_INCLUDE
        BLAS_LIBS=$BLAS_LIBS
        BLAS_LDFLAGS=$BLAS_LIBDIR
        BLAS_FOUND=yes],
       [AC_MSG_RESULT([no])
        AC_MSG_ERROR([no BLAS found])])
     LIBS=${blas_save_LIBS}
     CFLAGS=${blas_save_CFLAGS}
     AC_LANG_POP([C])

else # no --with-blas-* specified, so try pkg-config

  AC_CHECK_PROG([HAVE_PKG_CONFIG],[pkg-config],[yes])
  if test "x$HAVE_PKG_CONFIG" = xyes; then
     AC_MSG_CHECKING([for BLAS via pkg-config])
     cblas_module=`pkg-config --list-all | $AWK '$1=="cblas" {print $1}'`
     AS_IF([test "x$cblas_module" == "xcblas"],
           [
     AC_LANG_PUSH([C])
     blas_save_CFLAGS=${CFLAGS}
     CFLAGS=`pkg-config --cflags cblas`
     blas_save_LIBS=${LIBS}
     LIBS=`pkg-config --libs cblas`
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                       [[double x[3]; double y[3]; double alpha; cblas_daxpy(3, alpha, x, 1, y, 1);]])],
       [AC_MSG_RESULT([yes])
        AC_DEFINE([HAVE_BLAS],1,[Define if you have a (non-MKL) BLAS installation])
        BLAS_CPPFLAGS=`pkg-config --cflags cblas`
        BLAS_LIBS=`pkg-config --libs-only-l cblas`
        BLAS_LDFLAGS=`pkg-config --libs-only-L cblas`
        BLAS_FOUND=yes],
       [AC_MSG_RESULT([no])])
     LIBS=${blas_save_LIBS}
     CFLAGS=${blas_save_CFLAGS}
     AC_LANG_POP([C])
           ], [
      AC_MSG_RESULT([no])
])
  fi
fi
  if test "x$BLAS_FOUND" != "xyes"; then   # still not found? try last resort guess
    AC_MSG_CHECKING([for BLAS in default locations])
     AC_LANG_PUSH([C])
     blas_save_LIBS=${LIBS}
     LIBS=-lcblas
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                       [[double x[3]; double y[3]; double alpha; cblas_daxpy(3, alpha, x, 1, y, 1);]])],
       [AC_MSG_RESULT([yes])
        AC_DEFINE([HAVE_BLAS],1,[Define if you have a (non-MKL) BLAS installation])
        BLAS_CPPFLAGS=""
        BLAS_LIBS="-lcblas"
        BLAS_LDFLAGS=""],
       [AC_MSG_RESULT([no])
        AC_MSG_ERROR([No BLAS available])])
     LIBS=${blas_save_LIBS}
     CFLAGS=${blas_save_CFLAGS}
     AC_LANG_POP([C])
  fi
fi

AC_SUBST(BLAS_LIBS)
AC_SUBST(BLAS_CPPFLAGS)
AC_SUBST(BLAS_LDFLAGS)


AC_MSG_CHECKING([for Boost])
AC_LANG_PUSH([C++])
boost_save_CPP=$CPPLAGS
CPPFLAGS=$BOOST_CPPFLAGS
boost_save_LDFLAGS=$LDFLAGS
LDFLAGS=$BOOST_LDFLAGS
boost_save_LIBS=$LIBS
LIBS="-lboost_filesystem -lboost_system"
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <boost/filesystem.hpp> 
                                  using namespace boost::filesystem;]],
                   [[return file_size("/dev/null") > 0;]])],
  [AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])
                          AC_MSG_ERROR([no Boost available])])
LIBS=$boost_save_LIBS
LDFLAGS=$boost_save_LDFLAGS
CPPFLAGS=$boost_save_CPPFLAGS
AC_LANG_POP([C++])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_MSG_CHECKING([for cblas_daxpby])
  AC_LANG_PUSH([C++])
  mkl_save_CXXFLAGS=$CXXFLAGS
  CXXFLAGS=$BLAS_CPPFLAGS
  mkl_save_LIBS=$LIBS
  LIBS="$BLAS_LDFLAGS $BLAS_LIBS"

  if test "x$use_mkl" = xyes; then
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
