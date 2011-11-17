//
//
// $Header$
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#ifndef IEEE_FP_SETUP_H
#define IEEE_FP_SETUP_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_FENV_H
#include <fenv.h>

#define TRAPPED_FES (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)
#endif

#ifdef __SSE__
#include <xmmintrin.h>

#define MM_TRAPPED  (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW)
#endif

void ieeeFPsetup();

#endif
