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


void error(char *format, ...);
void coredump(char *format, ...);
void warning(char *format, ...);
void ensure(bool condition,char *errorIfFail, ...);


#endif
