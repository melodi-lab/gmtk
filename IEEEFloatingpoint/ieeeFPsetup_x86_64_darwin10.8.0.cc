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

#ifdef __SSE__
#include <string.h>
#include <xmmintrin.h>
#endif

#include <signal.h>
#include <sys/signal.h>
#include <fenv.h>

#ifdef __SSE__
#include "general.h"
#endif
#include "ieeeFPsetup_x86_64_darwin10.8.0.h"
#include "error.h"

#ifdef __SSE__
VCID(HGID)
#endif

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

struct sigaction act;

void ieeeFPsetup()
{
#ifdef __SSE__
	// printf("SSE\n");
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~MM_TRAPPED);
#else
	// printf("x87\n");
  feclearexcept(FE_ALL_EXCEPT);
  fenv_t env;
  if (fegetenv(&env)) {
    fprintf(stderr, "fgetenv() failed\n");
    abort();
  }
  env.__control &= ~(TRAPPED_FES);
  if (fesetenv(&env)) {
    fprintf(stderr,"fsetenv() failed\n");
    abort();
  }
#endif
  act.sa_sigaction = sighandler;
  act.sa_flags = SA_SIGINFO;
  if (sigaction(SIGFPE, &act, NULL) != 0) {
    fprintf(stderr,"Unable to install floating point signal handler\n");
    abort();
  }
}
