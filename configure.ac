#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([GMTK],m4_esyscmd([printf %s `head -n 1 RELEASE`]),
        [https://lungs.ee.washington.edu/trac/gmtk/newticket])
AM_INIT_AUTOMAKE([1.11 tar-ustar silent-rules -Wall -Werror foreign])

AC_ARG_ENABLE([gmtkViz], 
              [AS_HELP_STRING([--disable-gmtkViz],
                              [disable gmtkViz @<:@default: no@:>@])],
              [gmtkViz=${enableval}], [gmtkViz=yes])

# The following is for ticket 192. It's somewhat of a kludge.
# If the build cpp supports "-x assembler-with-cpp" (it's 
# likely only GCC's cpp does), we make GMTK's default cpp 
# command "$CPP -x assembler-with-cpp" in order to eliminate
# warnings about GMTK ASCII input files that aren't compatible 
# with the C lexical rules normally enforced by cpp - this 
# assumes the build cpp will be available and on the $PATH of 
# the end-user's machine. Otherwise, the default is just "cpp". 

# The user can override this at build-time with the 
# "--with-cppCmd" configure argument, or at GMTK run-time
# with the $GMTK_CPP_CMD environment variable.

# We're going to fiddle with CPP and its arguments, so
# remember what to reset them to when we're done
AC_REQUIRE_CPP
AC_LANG_PUSH([C])
gmtk_save_CPPFLAGS=${CPPFLAGS}
CPPFLAGS="-x assembler-with-cpp $CPPFLAGS"
#CPPFLAGS="-std=c89 $CPPFLAGS"

# see if the new flags handle the ticket 192 problems
AC_PREPROC_IFELSE([
  AC_LANG_SOURCE([[
#define CONCAT(a,b) a##b
#define CONCAT2(a,b) CONCAT(a,b)
CONCAT2(dt/,foo)
% GMTK's comments are not compatible with cpp]])],
  [gmtkCPP="$CPP -x assembler-with-cpp"],
  [gmtkCPP="cpp"])

# restore original environment
CPPFLAGS=$gmtk_save_CPPFLAGS
AC_LANG_POP([C])


AC_ARG_WITH([cppCmd],
            [AS_HELP_STRING([--with-cppCmd],
               [specify C pre-processor command GMTK uses to process ASCII input files @<:@default: "gcc -E -x assembler-with-cpp" if the build cpp is "gcc -E", otherwise just "cpp"@:>@])],
            [],
            [with_cppCmd=$gmtkCPP])

# end ticket 192 stuff

AC_ARG_ENABLE([gzip],
              [AS_HELP_STRING([--disable-gzip],
                              [disable piping .gz files through gzip @<:@default: no@:>@])],
              [gzip=${enableval}], [gzip=yes])
AC_ARG_ENABLE([bzip2],
              [AS_HELP_STRING([--disable-bzip2],
                              [disable piping .bz2 files through bzip2 @<:@default: no@:>@])],
              [bzip2=${enableval}], [bzip2=yes])
AC_ARG_ENABLE([cppIfASCII],
	      [AS_HELP_STRING([--disable-cppIfASCII],
                              [disable piping ASCII parameter files through cpp @<:@default: no@:>@])],
              [cppIfASCII=${enableval}], [cppIfASCII=yes])
AC_ARG_ENABLE([tmpLocalPool],
              [AS_HELP_STRING([--disable-tmpLocalPool],
                              [disable use of temporary local clique value pool @<:@default: no@:>@])],
              [tmpLocalPool=${enableval}], [tempLocalPool=yes])

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
              [the smallest log value that we represent (note that values > -1e16 may cause arithmetic problems (0+0 != 0) @<:@default: -1.0E17@:>@])],
	    [],
            [with_LZERO=(-1.0E17)])

AC_ARG_ENABLE([assertions],
              [AS_HELP_STRING([--disable-assertions],
                              [disable assertion checking @<:@default: no@:>@])],
              [assertions=${enableval}], [assertions=yes])

AC_ARG_ENABLE([model-checks],
              [AS_HELP_STRING([--disable-model-checks],
                              [disable run-time checks for models @<:@default: no@:>@])],
              [model_checks=${enableval}], [model_checks=yes])

AC_ARG_ENABLE([continuous-cardinality-warning],
              [AS_HELP_STRING([--disable-continous-cardinality-warning],
                              [disable warning about continous variables with non-zero cardinality in trifiles written by older versions of GMTK @<:@default: no@:>@])],
              [continuous_cardinaliy_warning=${enableval}], [continuous_cardinality_warning=yes])

# Checks for programs.

AC_PROG_AWK

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


# Autotest suite
AC_CONFIG_TESTDIR([tests],[tksrc:featureFileIO:miscSupport])
AC_CONFIG_FILES([tests/Makefile tests/atlocal])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([IEEEFloatingpoint
                   miscSupport
                   featureFileIO
                   mitk
                   tksrc])
AC_OUTPUT
