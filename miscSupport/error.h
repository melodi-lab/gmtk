/*
    $Header$
  
    Simple fatal error function.
    Jeff Bilmes <bilmes@cs.berkeley.edu>
*/


#ifndef ERROR_H
#define ERROR_H

#ifndef EXIT_SUCCESS
#include <stdlib.h>
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#define EXIT_FAILURE (1)
#endif
#endif

void error(const char * const format, ...);
void coredump(const char * constformat, ...);
void warning(const char * const format, ...);
void ensure(const bool condition,const char * const errorIfFail, ...);


#endif
