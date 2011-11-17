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


#if HAVE_CONFIG_H
#include <config.h>
#endif

#include<stdio.h>
#if HAVE_FENV_H
#include <fenv.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "hgstamp.h"
#include "general.h"
#include "ieeeFPsetup_x86_64_linux-gnu.h"
#include "error.h"

VCID(HGID)

void sighandler(int signum, siginfo_t *info, void *ptr) 
{
  if (signum == SIGFPE) {
    fprintf(stderr,"Received floating point (FP) exception: ");
    int flags = info->si_code;
    bool found = false;

    if (flags == FPE_INTDIV) {
      found = true;
      fprintf(stderr, "integer division by 0\n");
    }
    if (flags == FPE_INTOVF) {
      found = true;
      fprintf(stderr, "integer overflow\n");
    }
    if (flags == FPE_FLTDIV) {
      found = true;
      fprintf(stderr, "division by 0\n");
    }
    if (flags == FPE_FLTOVF) {
      found = true;
      fprintf(stderr, "overflow\n");
    }
    if (flags == FPE_FLTUND) {
      found = true;
      fprintf(stderr, "underflow\n");
    }
    if (flags == FPE_FLTRES) {
      found = true;
      fprintf(stderr, "inexact result\n");
    }
    if (flags == FPE_FLTINV) {
      found = true;
      fprintf(stderr, "invalid operation\n");
    }
    if (flags == FPE_FLTSUB) {
      found = true;
      fprintf(stderr, "subscript out of range\n");
    }
    if (!found) {
      fprintf(stderr, "Can't determine FP exception type %d\n", flags);
    }
    fprintf(stderr,"Process purposely exiting with a core dump due to FP Exception....\n");
    abort(); 
  } else {
    fprintf(stderr,"Caught signal %d, returning\n",signum);    
  }
}


void ieeeFPsetup()
{
#if HAVE_FENV_H
  feclearexcept(FE_ALL_EXCEPT); // clear all exception status bits
  feenableexcept(TRAPPED_FES);  // turn on trapping for our FEs
#endif
  struct sigaction act;
  memset(&act, 0, sizeof(act));
  act.sa_sigaction = sighandler;
  act.sa_flags = SA_SIGINFO;
  if (sigaction(SIGFPE, &act, NULL) != 0) {
    fprintf(stderr,"Unable to install floating point signal handler\n");
    abort();
  }
}
