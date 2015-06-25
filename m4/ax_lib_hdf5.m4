# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HDF5
#
# DESCRIPTION
#
#   This macro provides tests of the availability of HDF5 library.
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to the HDF5 helper script h5cc or h5pcc.
#
#   If HDF5 is successfully found, this macro calls
#
#     AC_SUBST(HDF5_VERSION)
#     AC_SUBST(HDF5_CXX)
#     AC_SUBST(HDF5_CXXFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)
#     AC_DEFINE(HAVE_HDF5)
#
#   and sets with_hdf5="yes".
#   H5CC will contain the appropriate wrapper script location.
#
#   If HDF5 is disabled or not found, this macros sets with_hdf5="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. HDF5_{CXX,CPP,LD}FLAGS may be used when building with C++.
#
#   To use the macro, one would code the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_HDF5()
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "HDF5 support:  $with_hdf5"
#
#   You could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that HDF5 uses:
#
#     AX_LIB_HDF5()
#     if test "$with_hdf5" = "yes"; then
#             CXX="$HDF5_CXX"
#     else
#             AC_MSG_ERROR([Unable to find HDF5, we need ${HDF5_CXX}.])
#     fi
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2015 Jeff Bilmes  (GMTK-specific modifications)
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 11

AC_DEFUN([AX_LIB_HDF5], [

AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(
    [--with-hdf5=[yes/no/PATH]], [location of h5c++ for HDF5 configuration]),
  [if test "$withval" = "no"; then
     with_hdf5="no"
   elif test "$withval" = "yes"; then
     with_hdf5="yes"
   else
     with_hdf5="yes"
     H5CC="$withval"
   fi],
   [with_hdf5="yes"]
)

dnl Set defaults to blank
HDF5_CXX=""
HDF5_VERSION=""
HDF5_CXXFLAGS=""
HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_LIBS=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_hdf5" = "yes"; then
    if test -z "$H5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [H5CC], [h5c++], [])
    else
        AC_MSG_CHECKING([Using provided HDF5 C++ wrapper])
        AC_MSG_RESULT([$H5CC])
    fi
    AC_MSG_CHECKING([for HDF5 libraries])
    if test ! -f "$H5CC" || test ! -x "$H5CC"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN(
[Unable to locate HDF5 compilation helper script 'h5c++'.
Please specify --with-hdf5=<LOCATION> as the full path to h5c++.
HDF5 support is being disabled (equivalent to --with-hdf5=no).])
        with_hdf5="no"
        with_hdf5_fortran="no"
    else
        dnl Get the h5cc output
        HDF5_SHOW=$(eval $H5CC -show)

        dnl Get the actual compiler used
        HDF5_CXX=$(eval $H5CC -show | $AWK '{print $[]1}')
        if test "$HDF5_CXX" = "ccache"; then
            HDF5_CXX=$(eval $H5CC -show | $AWK '{print $[]2}')
        fi

        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "HDF5 Version: X.Y.Z"
        HDF5_VERSION=$(eval $H5CC -showconfig | $GREP 'HDF5 Version:' \
            | $AWK '{print $[]3}')

        dnl A ideal situation would be where everything we needed was
        dnl in the AM_* variables. However most systems are not like this
        dnl and seem to have the values in the non-AM variables.
        dnl
        dnl We try the following to find the flags:
        dnl (1) Look for "NAME:" tags
        dnl (2) Look for "H5_NAME:" tags
        dnl (3) Look for "AM_NAME:" tags
        dnl
        HDF5_tmp_flags=$(eval $H5CC -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

        dnl Find the installation directory and append include/
        HDF5_tmp_inst=$(eval $H5CC -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

        dnl Add this to the CPPFLAGS
        HDF5_CPPFLAGS="-I${HDF5_tmp_inst}/include"

        dnl Now sort the flags out based upon their prefixes
        for arg in $HDF5_SHOW $HDF5_tmp_flags ; do
          case "$arg" in
            -I*) echo $HDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_CPPFLAGS="$arg $HDF5_CPPFLAGS"
              ;;
            -L*) echo $HDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LDFLAGS="$arg $HDF5_LDFLAGS"
              ;;
            -l*) echo $HDF5_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LIBS="$arg $HDF5_LIBS"
              ;;
          esac
        done

        HDF5_LIBS="$HDF5_LIBS -lhdf5_cpp -lhdf5"
        AC_MSG_RESULT([yes (version $[HDF5_VERSION])])

	AC_SUBST([HDF5_VERSION])
        AC_SUBST([HDF5_CXX])
	AC_SUBST([HDF5_CXXFLAGS])
	AC_SUBST([HDF5_CPPFLAGS])
	AC_SUBST([HDF5_LDFLAGS])
	AC_SUBST([HDF5_LIBS])

        AC_MSG_CHECKING([for HDF5 >= 1.8 C++ bindings])
        AC_LANG_PUSH([C++])
        hdf_check_save_CPP=${CPPFLAGS}
        CPPFLAGS="${HDF5_CPPFLAGS}"
        hdf_check_save_LIBS=${LIBS}
        LIBS="${HDF5_LIBS} -lhdf5_cpp"
        hdf_check_save_LD=${LDFLAGS}
        LDFLAGS="${HDF5_LDFLAGS}" 
        # H5Lregister should only be present in >= 1.8
        AC_LINK_IFELSE([AC_LANG_PROGRAM([
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
],[
  H5Lregister(NULL); 
  H5File file("/dev/null", H5F_ACC_RDONLY);
  hsize_t dims@<:@2@:>@; hsize_t maxdims@<:@2@:>@; 
  DataSpace *discDataspace = new DataSpace(2, dims, maxdims); 
  DSetCreatPropList prop;
  file.createDataSet(H5std_string("/dev/null"), PredType::NATIVE_UINT32, *discDataspace, prop)])],
           AC_SUBST([HAVE_HDF5],[yes]),[])
        LDFLAGS=${hdf_check_save_LDFLAGS}
        LIBS=${hdf_check_save_LIBS}
        CPPFLAGS=${hdf_check_save_CPPFLAGS}
        AC_LANG_POP([C++])
        if test x"$HAVE_HDF5" = x"yes"; then 
          CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
          LIBS="$LIBS $HDF5_LIBS -lhdf5_cpp"
          LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
          AC_DEFINE([HAVE_LIBHDF5_CPP],[1],[Are the C++ bindings for HDF5 available?])
          AC_MSG_RESULT([yes])
        else
          AC_MSG_RESULT([no])
          AC_MSG_WARN(
[HDF5 library installation was found but not usable. 
Be sure HDF5 version is >= 1.8, C++ bindings are installed, 
and use the same compiler to build HDF5 and GMTK. 
HDF5 support is being disabled (equivalent to --with-hdf5=no).])
        fi
	AC_DEFINE([HAVE_HDF5], [1], [Defined if you have HDF5 support])
    fi
fi
])
