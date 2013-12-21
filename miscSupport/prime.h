#ifndef PRIME_H
#define PRIME_H

/*  
 * Written by Richard Rogers <rprogers@ee.washington.edu> 
 *
 * Copyright (C) 2013 Richard Rogers 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_INTTYPES_H
   /* The ISO C99 standard specifies that the macros in inttypes.h must
      only be defined if explicitly requested. */
#  ifndef __STDC_FORMAT_MACROS
#    define __STDC_FORMAT_MACROS 1
#  endif
#  include <inttypes.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif

#ifndef true
#  define true 1
#endif
#ifndef false
#  define false 0
#endif

#ifdef __cplusplus
extern "C" {
#endif

int prime32(uint32_t n);

#ifdef __cplusplus
}
#endif

#endif
