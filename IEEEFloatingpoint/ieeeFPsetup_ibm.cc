//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// Solaris (and quite possibly very specific to a particular version
// of Solaris).
// 
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu



#include <stdio.h>
#include <string.h>


#include <signal.h>
#include <sys/signal.h>

#include <siginfo.h>
#include "general.h"
#include "ieeeFPsetup.h"
#include "error.h"


VCID("$Header$");

#ifdef HAVE_NONSTANDARD_ARITHMETIC

#endif

void fp_sigaction(int sig, siginfo_t* s_info, void *foo)
{
}

void ieeeFPsetup()
{

}





