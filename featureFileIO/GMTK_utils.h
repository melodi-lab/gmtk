//static char rcsid = "$Id";

#ifndef GMTK_UTILS_INCLUDED
#define GMTK_UTILS_INCLUDED

#include <stdio.h>
#include "error.h"

// read line from ascii file, check line length and empty lines

extern size_t fileSize(FILE *f);
extern size_t checklines(char *buf, char **pointers);
extern int emptyline(char *p1, char *p2);
extern char *skip_whitespace(char *);
extern char *skip_chars(char *);

#endif
