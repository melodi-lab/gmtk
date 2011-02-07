//
//
// Setup code so that FPE's print out not only that an FPE occured
// but what kind of FPE occured. This code is very specific to
// Solaris (and quite possibly very specific to a particular version
// of Solaris).
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

void fp_sigaction(int sig, siginfo_t* s_info, void *foo)
{
}

void ieeeFPsetup()
{

}





