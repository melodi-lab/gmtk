#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`), 
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([GMTK_DlinkMatrix.h])
AC_CONFIG_HEADERS([config.h])

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

AC_ARG_ENABLE([tmpLocalPool],
              [AS_HELP_STRING([--disable-tmpLocalPool],
                              [disable use of temporary local clique value pool @<:@default: no@:>@])],
              [tmpLocalPool=${enableval}], [tmpLocalPool=yes])
if test x"$tmpLocalPool" = x"yes"; then
  AC_DEFINE([USE_TEMPORARY_LOCAL_CLIQUE_VALUE_POOL],[1],[Hold clique values until after pruning])
fi

AC_ARG_ENABLE([continuous-cardinality-warning],
              [AS_HELP_STRING([--disable-continous-cardinality-warning],
                              [disable warning about continous variables with non-zero cardinality in trifiles written by older versions of GMTK @<:@default: no@:>@])],
              [continuous_cardinaliy_warning=${enableval}], [continuous_cardinality_warning=yes])
if test x"$continuous_cardinality_warning" = x"yes"; then
  AC_DEFINE([CONTINUOUS_CARDINALITY_WARNING],[1],[Warn about non-zero continuous RV cardinality in trifile])
fi

AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

AM_CONDITIONAL([NO_ASSERTIONS],[test x$assertions != xyes])
AM_CONDITIONAL([NO_MODEL_CHECKS],[test x$model_checks != xyes])

AC_ARG_WITH([LZERO],
            [AS_HELP_STRING([--with-LZERO],
              [the smallest log value that we represent (note that values > -1e16 may cause arithmetic problems (0+0 != 0) @<:@default: -1.0E17@:>@])],
	    [],
            [with_LZERO=(-1.0E17)])

badLZERO=`echo ${with_LZERO} | awk '{if ($1 > -1.0e16) {print "bad"}}'`
if test x"$badLZERO" = x"bad"; then
  AC_WARN([Values for LZERO greater than -1.0e16 may cause an arithmetic issue so that 0+0 != 0])
fi 
AC_DEFINE_UNQUOTED([LZERO], [(${with_LZERO})], [smallest log value that we represent])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LEX
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

AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])
if test x"$gmtkViz" = x"yes"; then
   AC_CHECK_PROG([HAVE_WXCONFIG],[wx-config],[yes])
   if test x"$HAVE_WXCONFIG" = x"yes"; then
     AC_LANG_PUSH([C++])
     wx_check_save_CXXFLAGS=${CXXFLAGS}
     CXXFLAGS=`wx-config --cxxflags`
     wx_check_save_LIBS=${LIBS}
     LIBS=`wx-config --libs`
     AC_LINK_IFELSE([AC_LANG_PROGRAM([#include "wx/wx.h"],
                                     [wxString s("foo");s.Capitalize()])],
        AC_SUBST([HAVE_WXWIDGETS],[yes]),
        AC_MSG_WARN(
	 [wx-config found but not wxWidgets >= 2.9. gmtkViz will not be built.])
     )
     LIBS=${wx_check_save_LIBS}
     CXXFLAGS=${wx_check_save_CXXFLAGS}
     AC_LANG_POP([C++])
   fi
fi
if test x"$HAVE_WXWIDGETS" = x"yes"; then
   AC_SUBST([GMTKVIZ],[gmtkViz])
   AC_SUBST([VIZCXXFLAGS],[`wx-config --cxxflags`])
   AC_SUBST([VIZLDFLAGS],[`wx-config --libs`])
fi

     AC_ARG_WITH(
       [ps2eps],
       [AS_HELP_STRING([--with-ps2eps],
          [specify default ps2eps command ('no' or '--without-ps2eps' generates PostScript) @<:@default: ps2eps@:>@])],
       [AS_IF([test "x$with_ps2eps" != "xno"],
          [AC_DEFINE_UNQUOTED([PS2EPS_CMD], ["${with_ps2eps}"],[default ps2eps command])
           AC_DEFINE([EPS_EXT], ["eps"], [extension for print to file])
           AC_DEFINE([PRINT2FILE_ABBREV], ["EPS"], [abbreviation for print to file format])
           AC_DEFINE([PRINT2FILE_FORMAT], ["Encapsulated PostScript"], [format for print to file])
          ],
          [AC_DEFINE([PS2EPS_CMD], ["cat"], [pass-thru raw PostScript])
           AC_DEFINE([EPS_EXT], ["ps"], [extension for print to file])
           AC_DEFINE([PRINT2FILE_ABBREV], ["PS"], [abbreviation for print to file format])
           AC_DEFINE([PRINT2FILE_FORMAT], ["PostScript"], [format for print to file])
          ])
       ],
       [AC_DEFINE([PS2EPS_CMD], ["ps2eps"], [default ps2eps command])
        AC_DEFINE([EPS_EXT], ["eps"], [extension for print to file])
        AC_DEFINE([PRINT2FILE_ABBREV], ["EPS"], [abbreviation for print to file format])
        AC_DEFINE([PRINT2FILE_FORMAT], ["Encapsulated PostScript"], [format for print to file])
       ]
     )

case "${host}" in
*apple*) 
  if test x"$HAVE_WXWIDGETS" = x"yes"; then
    AC_SUBST([OSXBUNDLE],[gmtkViz.app])

#   The EnableTop() bug (ticket 64) is fixed in wxWidgets 2.9.2
#   and gmtkViz souce checks if ENABLE_TOP_OFFSET is
#   needed using wxCHECK_VERSION(2,9,2)
    AC_DEFINE([GMTK_WX_OSX],[1],[Do we need to activate OSX fixes?])
  fi ;;
*-sun-solaris*) AC_CHECK_LIB([sunmath],[annuity]) ;;
*cygwin*) AC_SUBST([XOPEN],[-D__USE_XOPEN2K]) ;;
esac
AM_CONDITIONAL([CYGWIN_FLAGS],[test x"$XOPEN" = x"-D__USE_XOPEN2K"])
AM_CONDITIONAL([MAKE_BUNDLE],[test  x$OSXBUNDLE = xgmtkViz.app])
AC_CHECK_LIB([m], [log])


AC_LANG_PUSH([C++])
hdf_check_save_LIBS=${LIBS}
LIBS="${LIBS} -lhdf5_cpp"
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif],
[H5File file("/dev/null", H5F_ACC_RDONLY )])],
        AC_SUBST([HAVE_HDF5],[yes]),[])
LIBS=${hdf_check_save_LIBS}
AC_LANG_POP([C++])
if test x"$HAVE_HDF5" = x"yes"; then 
  LIBS="-lhdf5_cpp $LIBS"
  AC_DEFINE([HAVE_LIBHDF5_CPP],[1],[Are the C++ bindings for HDF5 available?])
fi


AC_ARG_WITH([detBrush],
            [AS_HELP_STRING([--with-detBrush],
              [select gmtkViz brush for deterministic nodes (none, cross, grey) @<:@default=none@:>@])],
            [],
            [with_detBrush=none])
case "${with_detBrush}" in
none)
  AC_DEFINE([NO_DET_BRUSH],[1],[no special brush for deterministic variables])
  ;;
cross)
  AC_DEFINE([CROSS_DET_BRUSH],[1],[cross-hatch brush for deterministic variables])
  ;;
grey)
  AC_DEFINE([GREY_DET_BRUSH],[1],[solid grey brush for deterministic variables])
  ;;
*)
  AC_MSG_ERROR([Must specify --with-detBrush= one of none, cross, grey]) 
  ;;
esac

dnl The dlopen() function is in the C library for *BSD and in
dnl libdl on GLIBC-based systems
AC_SEARCH_LIBS([dlopen], [dl dld], [], [
  AC_MSG_ERROR([unable to find the dlopen() function])
])

# Checks for header files.
AC_PATH_XTRA
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h dlfcn.h])

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
AC_CHECK_FUNCS([alarm floor memset pow regcomp sqrt strchr strerror strstr strtol log1p getline])
AC_FUNC_FSEEKO

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
