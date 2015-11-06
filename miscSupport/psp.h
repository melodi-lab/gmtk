/*-
 * psp.h
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

#ifndef PSP_H
#define PSP_H

#include <stdio.h>

void psp2(FILE*f,const int numSpaceChars,const char c = ' ');

void psp(FILE*f,const int numSpaceChars,const char c = ' ');

#endif
