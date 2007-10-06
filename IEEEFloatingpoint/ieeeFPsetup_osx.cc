//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// OSX (and quite possibly very specific to a particular version
// of OSX).
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


void fp_sigaction(int sig, siginfo_t* s_info, void *foo)
{
}

void ieeeFPsetup()
{

}





