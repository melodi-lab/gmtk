#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# autoconf initialization
AC_PREREQ(2.59)
AC_INIT(GMTK, 0.0, gmtk-bugs@ssli.ee.washington.edu)
AC_CONFIG_AUX_DIR(config)

# This identifies the target CPU/OS to select the source files containing
# the correct platform-specific FPU code for the target. The source files
# are named ieeeFPsetup_$(IEEEFLAVOR).{h,cc} in the Makefile.

# TODO: rename IEEE support files ieeeFPsetup_$(target_cpu)_$(target_os).{h,cc}
#       or maybe $(target_fpu)_$(target_os) [would have to invent target_fpu]
#       The code really depends on both OS and (C|F)PU, but x86_64 and i*86 probably
#       both use the same x87 FPU code. I wonder if x86_64 CPUs show up as
#       i*86 $(target_cpu) when running 32-bit OSes?
AC_CANONICAL_TARGET
case "$target" in
  x86_64*-linux* | i*86*-linux*)
    IEEEFLAVOR=linux
    ;;
  x86_64*-darwin* | i*86*darwin*)
    IEEEFLAVOR=osx
    ;;
  *-cygwin*)
    IEEEFLAVOR=cygwin
    ;;
  *)
    IEEEFLAVOR=${target_cpu}_${target_os}
    # TODO: could create a make target to copy a do-nothing default implementation to 
    # the appropriate filenames, provide instructions to implement the files
    AC_MSG_WARN([Create files ieeeFPsupport_${target_cpu}_${target_os}.{h,cc} with code for the target platform's FPU.])
    AC_MSG_ERROR([Unsupported CPU/OS ${target_cpu}/${target_os}.])
    ;;
esac
AC_SUBST(IEEEFLAVOR) 
AM_INIT_AUTOMAKE

# Check that we're in the correct directory by verifying the existence of the file
AC_CONFIG_SRCDIR([ieeeFPsetup_linux.cc])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
