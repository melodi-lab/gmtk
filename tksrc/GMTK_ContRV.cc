/*
 * GMTK_ContRV.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (C) 2001 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */



#include "general.h"
#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)


#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_ContRV.h"
#include "debug.h"


/*-
 *-----------------------------------------------------------------------
 * printSelf()
 *      prints a one-line summary of the detailed information about this RV.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
void ContRV::printSelf(FILE *f,bool nl)
{
  printNameFrameValue(f,false);
  fprintf(f," continuous%s",nls(nl));
}



/*-
 *-----------------------------------------------------------------------
 * cloneRVShell()
 *      clones a shell of the current random variable (see GMTK_RV.h for docs)
 *
 * Preconditions:
 *      RV must be instantiated and with parameters (i.e., what lives in the template RVs).
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      self is printed.
 *
 *-----------------------------------------------------------------------
 */
ContRV* ContRV::cloneRVShell()
{
  ContRV*rv = (ContRV*)RV::cloneRVShell();
  rv->conditionalMixtures = conditionalMixtures;
  rv->curMappingOrDirect = &rv->conditionalMixtures[0];
  return rv;
}


