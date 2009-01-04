/*
    $Header$
  
    Simple fatal error function.
    Jeff Bilmes <bilmes@cs.berkeley.edu>
    $Header$
*/


#ifndef __GNUC__
enum bool { false = 0, true = 1 };
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "error.h"


void
error(const char * const format, ...)
{
  va_list ap;
  va_start(ap,format);
  /* print out remainder of message */
  (void) vfprintf(stderr, format, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
#ifdef ERROR_DOES_ABORT
  (void) abort();
#else
  (void) exit(EXIT_FAILURE);
#endif
}

void
coredump(const char * const format, ...)
{
  va_list ap;
  va_start(ap,format);
  /* print out remainder of message */
  (void) vfprintf(stderr, format, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
  (void) abort();
}

void
warning(const char * const format, ...)
{
  va_list ap;
  va_start(ap,format);
  /* print out remainder of message */
  (void) vfprintf(stderr, format, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
}

void
ensure(const bool condition,const char * const errorIfFail, ...)
{
  if (!condition) {
    va_list ap;
    va_start(ap,errorIfFail);
    /* print out remainder of message */
    (void) vfprintf(stderr, errorIfFail, ap);
    va_end(ap);
    (void) fprintf(stderr, "\n");
    (void) exit(EXIT_FAILURE);
  }
}


#ifdef MAIN

int main()
{
  warning("This is a warning with output %d %f (%s)\n",
	  4,4.5,"A string");
  error("This is a fatal error with output %d %f (%s), program should die after this.\n",
	  4,4.5,"A string");
  return 0;
}

#endif
