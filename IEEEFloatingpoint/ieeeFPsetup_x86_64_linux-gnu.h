/*
**
** This code will set up the FPU on x86 archiectures to trap with
** SIGFPEs when an inf or nan occurs (or some other fp exception depending
** on what is uncommented below).
**
** The code is very machine specific to Linux (and probably to a particular
** version of Linux).
**
** Jeff Bilmes
** <bilmes@ee.washington.edu> 
**
**
*/

#ifndef IEEESETUP_H
#define IEEESETUP_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_FENV_H
#include <fenv.h>

// enable trapping on the following FPU exceptions
#define TRAPPED_FES (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)
#endif

extern void ieeeFPsetup();

#endif

