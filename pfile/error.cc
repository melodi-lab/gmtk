/*
    $Header$
  
    Simple fatal error function.
    Jeff Bilmes <bilmes@cs.berkeley.edu>
    $Header$
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


#include "error.h"

void
error(char *format, ...)
{
  va_list ap;
  va_start(ap,format);
  /* print out remainder of message */
  (void) vfprintf(stderr, format, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
  (void) exit(EXIT_FAILURE);
}
