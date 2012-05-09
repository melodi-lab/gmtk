//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// CPU and OS.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include <stdio.h>
#include <string.h>


#include <signal.h>
#include <sys/signal.h>
#include "general.h"
#ifdef warning
#warning "Change the ANY_ANY below to $(host_cpu)_$(host_os)" 
#endif
#include "ieeeFPsetup_ANY_ANY.h"
#include "error.h"


VCID(HGID)


void ieeeFPsetup()
{
  // does not yet do anything for this architecture.
}





