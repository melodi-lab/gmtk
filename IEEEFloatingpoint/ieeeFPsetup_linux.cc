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



#include<stdio.h>
#include<fpu_control.h>
#include <signal.h>
#include <stdlib.h>

#include "ieeeFPsetup.h"


/*
** Gets the i386 status bits
*/
#define _FPU_GETSTW(sw) __asm__ ("fnstsw %0" : "=m" (*&sw))

void sighandler(int sigarg) 
{
  // unsigned int cw=0;
  if (sigarg == SIGFPE) {
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
    fprintf(stderr,"Caught signal %d, returning\n",sigarg);    
  }
}


void ieeeFPsetup()
{
  unsigned int cw=0;
  // unsigned int sw=0;

  /* set the signal handler */
  signal(SIGFPE,(void(*)(int))sighandler);


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

}




