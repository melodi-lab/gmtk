#ifndef POPEN_ERR_H
#define POPEN_ERR_H

/*
 * popen_err.h
 * 
 * A version of popen() that also makes stderr available
 * so that warnings/error messages can be prefixed.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 * 
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define POPEN_BUF_SIZE 4095
#define POPEN_MAX_ARGC 1023

/*
 * Normal popen() only gives us access to the forked process' stdout.
 * To address https://j.ee.washington.edu/trac/gmtk/ticket/389 we need
 * the process' stderr as well. So popen_err() gives us access to both.
 * Any lines of output on the process' stderr are prefixed by prefix.
 */

FILE *popen_err(char const *command, char const *type, char const *prefix);

int pclose_err(FILE *stream);

#ifdef __cplusplus
}
#endif

#endif
