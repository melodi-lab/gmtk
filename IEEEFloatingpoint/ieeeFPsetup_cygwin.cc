//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// Cygwin.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu


#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>


#include <signal.h>
#include <sys/signal.h>

#include "hgstamp.h"
#include "general.h"
#include "ieeeFPsetup.h"
#include "error.h"


VCID(HGID)


void ieeeFPsetup()
{
  // does not yet do anything for this architecture.
}





