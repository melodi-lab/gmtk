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


// Explicitly define '__EXTENSIONS__' so that struct sigaction is defined
// in the following include files.
#ifndef __EXTENSIONS__
#define __EXTENSIONS__
#endif

#include <signal.h>
#include <sys/signal.h>

#include <siginfo.h>
#include "general.h"
#include "ieeeFPsetup.h"
#include "error.h"


VCID("$Header$");

#ifdef HAVE_NONSTANDARD_ARITHMETIC
extern "C" void nonstandard_arithmetic();
extern "C" void ieee_retrospective(FILE*);
#endif

void fp_sigaction(int sig, siginfo_t* s_info, void *foo)
{
  char *errstr;
  switch (s_info->si_code) {
  case FPE_INTDIV:      
    errstr = "integer divide by zero"; break;
  case FPE_INTOVF:
    errstr = "integer overflow"; break;
  case FPE_FLTDIV:
    errstr = "floating point divide by zero"; break;
  case FPE_FLTOVF:
    errstr = "floating point overflow"; break;
  case FPE_FLTUND:
    errstr = "floating point underflow"; break;
  case FPE_FLTRES:
    errstr = "floating point inexact result"; break;
  case FPE_FLTINV:
    errstr = "invalid floating point operation"; break;
  case FPE_FLTSUB:
    errstr = "subscript out of range"; break;
  default:
    errstr = "unknown code"; break;
  }
  coredump("Program received signal SIGFPE (%d), Arithmetic exception (%s)",sig,errstr);
}

void ieeeFPsetup()
{
  
  struct sigaction act;
  act.sa_sigaction = fp_sigaction;
  ::memset(&(act.sa_mask),0,sizeof(act.sa_mask));
  act.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE,&act,NULL);
  fpsetmask(
	    FP_X_INV        /* invalid operation exception */
	    | FP_X_OFL      /* overflow exception */
	    // | FP_X_UFL      /* underflow exception */
	    | FP_X_DZ       /* divide-by-zero exception */
	    // | FP_X_IMP      /* imprecise (loss of precision) */
	    );

#ifdef HAVE_NONSTANDARD_ARITHMETIC
  // This presumably sets a bit in the FPU that keeps
  // denormals from traping and being handled in sofware.
  // Instead, (again presumably), denormals are truncated to zero.
  // WARNING: This should probably not be used when debugging.
  // On the ohter hand, this can significantly speed up a program 
  // w/o changing the results much.
  // Note: you need to have the sunmath library (-lsunmat) to
  // have this routine. This is part of the SUNWspro SC sun
  // optimizing compiler.
  nonstandard_arithmetic();

  // print out a status of the FPU bits.
  // ieee_retrospective(stdout); 
#endif

  // these are the options. Default is FP_RN
  // FP_RN         /* round to nearest representative number */
  // FP_RP         /* round to plus infinity */
  // FP_RM         /* round to minus infinity */
  // FP_RZ         /* round to zero (truncate) */
  // fpsetround(FP_RZ);

}





