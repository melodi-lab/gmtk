//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// Cygwin.
// 
// Written by: Jeff Bilmes
//             bilmes@ee.washington.edu



#include <stdio.h>
#include <string.h>


#include <signal.h>
#include <sys/signal.h>

#include "general.h"
#include "ieeeFPsetup.h"
#include "error.h"


VCID("$Header$")


void ieeeFPsetup()
{
  // does not yet do anything for this architecture.
}





