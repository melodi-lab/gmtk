/*-
 * GMTK_Signals.cc
 *     Handling of Unix signals
 *
 * Written by Chris Bartels & Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

#include "GMTK_Signals.h" 
#include "error.h"
#include "general.h"

#include <stdio.h>     /* standard I/O functions                         */
#include <signal.h>    /* signal name macros, and the signal() prototype */

// can't call this 'terminate' because there is a terminate() function
// in the std namespace
volatile bool sigterminate = false;


/*-
 *-----------------------------------------------------------------------
 * catch_sigusr1 
 *   The signal handler for SIGUSR1
 *
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   The sigterminate flag is set
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void catch_sigusr1(int sig_num)
{
  // Display warning 
  fprintf(stderr, "GMTK received sigusr1, terminating...\n");
  sigterminate = true;
}


/*-
 *-----------------------------------------------------------------------
 * catch_sigxcpu
 *   The signal handler for SIGXCPU
 *
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   program exits with status EXIT_RESOURCES_EXCEEDED
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void catch_sigxcpu(int sig_num)
{
  // Display warning 
  fprintf(stderr, "GMTK received SIGXCPU (maximum CPU time has been exceeded), terminating...\n");
  exit_program_with_status(EXIT_RESOURCES_EXCEEDED);

  //  this method requires testing TerminateSignalReceived() inside loops
  // sigterminate = true;
  //  in future, could use this method to allow cleaner exit, or
  //  possibly to limit the CPU time *per utterance* processed

}



/*-
 *-----------------------------------------------------------------------
 * InstallSignalHandlers 
 *   Installs signal handlers (currently only SIGUSR1 is caught).
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   Signal handlers are installed, sigterminate flag set to false
 *
 * Side Effects:
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void InstallSignalHandlers()
{
  sigterminate = false;

  signal(SIGUSR1, catch_sigusr1);
  signal(SIGXCPU, catch_sigxcpu);
}

/*-
 *-----------------------------------------------------------------------
 * TerminateSignalReceived 
 *   Tells if a signal has been received that indicates that the program 
 *   should terminate.  
 *
 * Preconditions:
 *   none
 * 
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   true if the program should terminate, false if not 
 *
 *-----------------------------------------------------------------------
 */
bool TerminateSignalReceived()
{
  return(sigterminate); 
}

