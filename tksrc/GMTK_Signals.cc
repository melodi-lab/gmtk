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

#include <stdio.h>     /* standard I/O functions                         */
#include <signal.h>    /* signal name macros, and the signal() prototype */

volatile bool terminate = false;

/*-
 *-----------------------------------------------------------------------
 * catch_sigusr1 
 *   The signal handler for SIGUSR1
 *
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   The terminate flag is set
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
  terminate = true;
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
 *   Signal handlers are installed, terminate flag set to false
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
  terminate = false;

  signal(SIGUSR1, catch_sigusr1);
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
  return(terminate); 
}

