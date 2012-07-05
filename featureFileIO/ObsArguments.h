/*
 * ObsArguments.h
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 *
 *   Defines most of the arguments for all obs- programs in one place. Defines all three of:
 *       1) the (static) variable names, if not defined elsewhere in GMTK.
 *       2) the argument names, definition, and documentation .
 *       3) code to check that the values of the arguments as specified on the command line is correct.
 * 
 *   The user of this file includes it three times, each time with one of the below defined.
 *       GMTK_ARGUMENTS_DEFINITION     // defines the argument as a static variable if needed
 *       GMTK_ARGUMENTS_DOCUMENTATION  // specifies documentation to arguments.h
 *       GMTK_ARGUMENTS_CHECK_ARGS     // emits code to check that the argument was set ok according to command line, and dies with error if not.
 * 
 *  Each argument is obtained by defining an appropriate #define with a number (0,1,2) with the arguments priority for arguments.h
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>

#ifdef HAVE_HG_H
#include "hgstamp.h"
#endif
#endif

#include "GMTK_Filter.h"

/*-----------------------------------------------------------------------------------------------------------*/
/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

/* initial definitions commonto all arguments */

#if defined(GMTK_ARGUMENTS_DEFINITION)

#elif defined(GMTK_ARGUMENTS_DOCUMENTATION)

#elif defined(GMTK_ARGUMENTS_CHECK_ARGS)

const char*const argerr = "ARG ERROR";
// include nop statement to avoid warning message.
(void)argerr;


#else
#endif


#include "GMTK_ObservationArguments.h"
