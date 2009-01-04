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

#include <xmmintrin.h>


VCID("$Header$")


void sighandler(int sigarg) 
{
  if (sigarg == SIGFPE) {
    unsigned int sw=0;
    int found=0;
    /* Ideally, this should work but it looks like for now
       the osx kernal is clearing the FPU status bits so we
       are not able, at this point, to figure out what type
       of floating point exception occured. Once the kernel
       is fixed, this code should start working.
    */
    sw = _MM_GET_EXCEPTION_STATE();

    fprintf(stderr,"Received floating point (FP) exception: ");
    if (sw & _MM_EXCEPT_INVALID) {
      fprintf(stderr,"invald operation\n");
      found = 1;
    }
    if (sw & _MM_EXCEPT_DENORM) {
      fprintf(stderr,"denormalized operand\n");
      found = 1;
    }
    if (sw & _MM_EXCEPT_DIV_ZERO) {
      fprintf(stderr,"divide by zero\n");
      found = 1;
    }
    if (sw & _MM_EXCEPT_OVERFLOW) {
      fprintf(stderr,"overflow\n");
      found = 1;
    }
    if (sw & _MM_EXCEPT_UNDERFLOW) {
      fprintf(stderr,"underflow\n");
      found = 1;
    }
    if (sw & _MM_EXCEPT_INEXACT) {
      fprintf(stderr,"(inexact) low precision\n");
      found = 1;
    }
    if (!found) {
      fprintf(stderr,"Can't determine FP exception type 0x%X\n",sw);
    }
    fprintf(stderr,"Process purposely exiting with a core dump due to FP Exception....\n");
    abort();
  } else {
    fprintf(stderr,"Caught signal %d, returning\n",sigarg);    
  }
}


void ieeeFPsetup()
{

  /* set the signal handler */
  signal(SIGFPE,(void(*)(int))sighandler);


  // a simple way of setting one mask
  // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

  unsigned int cw = _MM_GET_EXCEPTION_MASK();

  // set the ones we want
  cw = cw & ~( 
	      _MM_MASK_INVALID            /* Invalid operation, i.e., NaNs */	    
	      /* | _MM_MASK_DENORM */	  /* Denormalized operand */		    
	      | _MM_MASK_DIV_ZERO	  /* Zero-divide */			    
	      | _MM_MASK_OVERFLOW	  /* Overflow  */			    
	      /* | _MM_MASK_UNDERFLOW */  /* Underflow  */			    
	      /* | _MM_MASK_INEXACT */    /* Precision (inexact result) */        
	       );
	      
  _MM_SET_EXCEPTION_MASK( cw );

}



/*
 * various bits of information for MacOS

Here's another pointer, from Shoaib Kamil:

  Perhaps this is a good answer:

http://stackoverflow.com/questions/247053/enabling-floating-point-interrupts-on-mac-os-x-intel

  The short version is: SSE code doesn't seem to honor exception
  handling and therefore the compiler must emit x87 code only; the flags
  need to be set using inline assembly because of the missing
  feenableexcept() on MacOS.

================================================================================	

On Linux, feenableexcept and fedisableexcept can be used to control
the generation of SIGFPE interrupts on floating point exceptions. How
can I do this on Mac OS X Intel?

Inline assembly for enabling floating point interrupts is provided in 

http://developer.apple.com/documentation/Performance/Conceptual/Mac_OSX_Numerics/Mac_OSX_Numerics.pdf

pp. 7-15, but only for PowerPC assembly.

Geoffrey Irving

--------------------------------------------------------------------------------	

On Mac OS X this is moderately complicated. OS X uses the SSE unit for
all FP math by default, not the x87 FP unit. The SSE unit does not
honor the interrupt options, so that means that in addition to
enabling interrupts, you need to make sure to compile all your code
not to use SSE math.

You can disable the math by adding "-mno-sse -mno-sse2 -mno-sse3" to
your CFLAGS. Once you do that you can use some inline assembly to
configure your FP exceptions, with basically the same flags as Linux.

short fpflags = 0x1332 // Default FP flags, change this however you want. 
asm("fnclex");
asm("fldcw _fpflags");

The one catch you may find is that since OS X is built entirely using
sse there may be uncaught bugs. I know there used to be a big with the
signal handler not passing back the proper codes, but that was a few
years ago, hopefully it is fixed now.  link|offensive?


Louis Gerbarg

	
================================================================================		

Exceptions for sse can be enabled using _MM_SET_EXCEPTION_MASK from
xmmintrin.h. For example, to enable invalid (nan) exceptions, do

#include <xmmintrin.h>
...
_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);



*/	


