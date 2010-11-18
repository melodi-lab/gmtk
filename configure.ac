#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK],m4_esyscmd([printf %s `head -n 1 RELEASE`]),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])

# Checks for programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([IEEEFloatingpoint
                   miscSupport
                   featureFileIO
                   mitk
                   tksrc])
AC_OUTPUT
