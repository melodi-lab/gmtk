/*-
 * psp.cc
 *
 *     print indentation
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#  include "hgstamp.h"
#endif
#include "general.h"

#include "psp.h"

VCID(HGID)

void
psp2(FILE*f, const int numSpaceChars, const char c) {
  int tmp = numSpaceChars;
  // prints a '|' every 4 chars. 
  if (!tmp) return;
  fprintf(f,"|"); tmp--;
  while (tmp > 4) {
    fprintf(f,"   |");
    tmp -= 4;
  }
  while (tmp--)
    fprintf(f,"%c",c);
}

void
psp(FILE*f, const int numSpaceChars, const char c) {
  int tmp = numSpaceChars;
  while (tmp--)
    fprintf(f,"%c",c);
}
