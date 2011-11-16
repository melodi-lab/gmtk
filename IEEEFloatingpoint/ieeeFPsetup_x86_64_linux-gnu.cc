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
#include <fenv.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "hgstamp.h"
#include "general.h"
#include "ieeeFPsetup_x86_64_linux-gnu.h"
#include "error.h"

VCID(HGID)

/*
** Gets the i386 status bits
*/
#define _FPU_GETSTW(sw) __asm__ ("fnstsw %0" : "=m" (*&sw))


void sighandler(int signum, siginfo_t *info, void *ptr) 
{
#if 1
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
#else
  // unsigned int cw=0;
  if (signum == SIGFPE) {
    unsigned int sw=0;
    int found=0;
    /* Ideally, this should work but it looks like for now
       the linux kernal is clearing the FPU status bits so we
       are not able, at this point, to figure out what type
       of floating point exception occured. Once the kernel
       is fixed, this code should start working.
    */
    _FPU_GETSTW(sw);

    fprintf(stderr,"Received floating point (FP) exception: ");
    if (sw & 0x1) {
      fprintf(stderr,"(inexact) low precision\n");
      found = 1;
    }
    if (sw & 0x2) {
      fprintf(stderr,"underflow\n");
      found = 1;
    }
    if (sw & 0x4) {
      fprintf(stderr,"overflow\n");
      found = 1;
    }
    if (sw & 0x8) {
      fprintf(stderr,"divide by zero\n");
      found = 1;
    }
    if (sw & 0x10) {
      fprintf(stderr,"denormalized operand\n");
      found = 1;
    }
    if (sw & 0x20) {
      fprintf(stderr,"invald operation\n");
      found = 1;
    }
    if (sw & 0x40) {
      fprintf(stderr,"ES (floating-point exception summary)\n");
      found = 1;
    }
    if (!found) {
      fprintf(stderr,"Can't determine FP exception type\n");
    }
    fprintf(stderr,"Process purposely exiting with a core dump due to FP Exception....\n");
    abort();
  } else {
    fprintf(stderr,"Caught signal %d, returning\n",signum);    
  }
#endif
}

struct sigaction act;

#define TRAPPED_FES FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW

void ieeeFPsetup()
{
#if 1
  feclearexcept(TRAPPED_FES);
  feenableexcept(TRAPPED_FES);

  memset(&act, 0, sizeof(act));
  act.sa_sigaction = sighandler;
  act.sa_flags = SA_SIGINFO;
  if (sigaction(SIGFPE, &act, NULL) != 0) {
    fprintf(stderr,"Unable to install floating point signal handler\n");
    abort();
  }
#else
  unsigned int cw=0;
  // unsigned int sw=0;

  /* set the signal handler */
  if (signal(SIGFPE,(void(*)(int))sighandler) == SIG_ERR) {
    fprintf(stderr,"Unable to install floating point signal handler\n");
    abort();
  }

  //  get the current FPU control word
  // _FPU_GETCW(cw);
  // _FPU_GETSTW(sw);

  /* printf("Before setting cw = 0x%X, sw = 0x%X\n",cw,sw); */

  _FPU_GETCW(cw);
  /* change the control word to catch these FP exceptions */
  cw = cw & ~( 
	      _FPU_MASK_IM          /* Invalid operation, i.e., NaNs */
	      /* | _FPU_MASK_DM */  /* Denormalized operand */
	      | _FPU_MASK_ZM        /* Zero-divide */
	      | _FPU_MASK_OM        /* Overflow  */
	      /* | _FPU_MASK_UM */  /* Underflow  */
	      /* | _FPU_MASK_PM */  /* Precision (inexact result) */
	      );

  /* set the control word now */
  _FPU_SETCW(cw);

  /*
  _FPU_GETSTW(sw);
  _FPU_GETCW(cw);
  printf("After setting cw = 0x%X, sw = 0x%X\n",cw,sw);
   printf("%f\n",foo(a,b));
  _FPU_GETCW(cw);
  _FPU_GETSTW(sw);
  */
#endif
}




