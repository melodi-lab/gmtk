#                                               -*- Autoconf -*-

#
# Written by Richard Rogers <rprogers@uw.edu>
#
# Copyright (C) 2015 Jeff Bilmes
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
#

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([GMTK], m4_esyscmd(printf %s `head -n 1 ../RELEASE`),
        [https://j.ee.washington.edu/trac/gmtk/newticket])

AC_CONFIG_SRCDIR([gmtk.h2m])

AM_INIT_AUTOMAKE([1.12 silent-rules -Wall -Werror foreign])
AM_SILENT_RULES([yes])

# Checks for programs.

AM_MISSING_PROG(HELP2MAN, help2man, $missing_dir)
AC_PROG_SED

#AC_PROG_INSTALL

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
