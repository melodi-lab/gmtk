#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK],m4_esyscmd([printf %s `head -n 1 RELEASE`]),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
AM_INIT_AUTOMAKE([1.11 silent-rules -Wall -Werror foreign])

AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])
AC_ARG_WITH([logp], 
            [AS_HELP_STRING([--with-logp],
              [select logp implementation (log, table, smart, renorm) @<:@default: log@:>@])],
            [],
            [with_logp=log])
AC_ARG_WITH([detBrush],
            [AS_HELP_STRING([--with-detBrush],
              [select gmtkViz brush for deterministic nodes (none, cross, grey) @<:@default: none@:>@])],
            [],
            [with_detBrush=none])

AC_ARG_WITH([LZERO],
            [AS_HELP_STRING([--with-LZERO],
              [the smallest log value that we represent @<:@default: -1.0E10@:>@])],
	    [],
            [with_LZERO=(-1.0E10)])

AC_ARG_ENABLE([assertions],
              [AS_HELP_STRING([--disable-assertions],
                              [disable assertion checking @<:@default: no@:>@])],
              [assertions=${enableval}], [assertions=yes])

AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

# Checks for programs.

AC_PROG_AWK

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


# Autotest suite
AC_CONFIG_TESTDIR([tests],[tksrc])
AC_CONFIG_FILES([tests/Makefile tests/atlocal])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([IEEEFloatingpoint
                   miscSupport
                   featureFileIO
                   mitk
                   tksrc])
AC_OUTPUT
