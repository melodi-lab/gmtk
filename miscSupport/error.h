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


#define E_INFO_INCREMENT            (10)

#define E_INFO_NOTHING              (0*E_INFO_INCREMENT)
#define E_INFO_E_INFO_LOW_PRI       (1*E_INFO_INCREMENT)
#define E_INFO_E_INFO_MED_PRI       (2*E_INFO_INCREMENT)
#define E_INFO_E_INFO_HIGH_PRI      (3*E_INFO_INCREMENT)
#define E_INFO_WARNING_LOW_PRI      (4*E_INFO_INCREMENT)
#define E_INFO_WARNING_MED_PRI      (5*E_INFO_INCREMENT)
#define E_INFO_WARNING_HIGH_PRI     (6*E_INFO_INCREMENT)

extern int e_info_level;

void error(char *format, ...);
void coredump(char *format, ...);
void warning(char *format, ...);
void ensure(bool condition,char *errorIfFail, ...);

inline void info(int fo,char*format, ...) {
  if (fo >= e_info_level) {
    va_list ap;
    va_start(ap,format);
    /* print out remainder of message */
    (void) vfprintf(stdout, format, ap);
    va_end(ap);
    (void) fprintf(stdout, "\n");
  }
}



#endif
