/*
 * GMTK_Sw_HidDiscRV.cc
 *
 * Observed discrete random variable.
 * 
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *  $Header$
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 * The top level GMTK random variable object for the RV class hierarchy.
 *
 *
 *
 */



#include "general.h"
VCID("$Header$");

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <string.h>

#include "GMTK_Sw_HidDiscRV.h"



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
void Sw_HidDiscRV::printSelf(FILE *f,bool nl)
{
  printNameFrameValue(f,false);
  fprintf(f,"switching hidden discrete cardinality = %d%s",cardinality,nls(nl));
}



/*-
 *-----------------------------------------------------------------------
 * printSelfVerbose()
 *      prints a multi-line verbose description of this RV.
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
void Sw_HidDiscRV::printSelfVerbose(FILE *f)
{
  fprintf(f,"Switching Hidden Discrete Random variable:\n");
  printNameFrameValue(f,true);
  fprintf(f,"From line %d in file %s\n",rv_info.fileLineNumber,rv_info.rvFileName.c_str());
  fprintf(f,"RV has cardinality = %d\n",cardinality);
}

#if 0

/*-
 *-----------------------------------------------------------------------
 * identicalStructureWith.
 *      Returns true if this rv has identical structure with that of other.
 *      "identical structure" means that the r.v. have the same
 *      number, type, and cardinality parents. If this returns
 *      true, then it will be valid to tie parameters between
 *      these two random variables. Note that this routine
 *      might need to change for each subclass of this class.
 * 
 * Preconditions:
 *      Both rvs must have parents filled in.
 *
 * Postconditions:
 *      If function returns true, then the variables have
 *      identical structure, otherwise not.
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      indicating boolean.
 *
 *-----------------------------------------------------------------------
 */
bool
Sw_HidDiscRV::identicalStructureWith(RV& other)
{

  if (!other.discrete())
    return false;
  if (!other.switching())
    return false;

  return SwDiscRV::identicalStructureWith(*(SwDiscRV*)&other);
}



/*-
 *-----------------------------------------------------------------------
 * tieParametersWith()
 *      Ties the parameters of 'this' with whatever those of 'other' are. 
 *      'other' and 'this' must be identical structuraly, if the
 *       'checkStructure' option is true.
 * 
 * Preconditions:
 *      other must be a fully instantiated RV with parameters, and 'this'
 *      and 'other' must be structurally identical (if arg is true)
 *
 * Postconditions:
 *      'this' has the identical _tied_ parameters with 'other'
 *
 * Side Effects:
 *      Changes the internal parameter data structures, but does not delete anything.
 *
 * Results:
 *      returns nothing.
 *
 *-----------------------------------------------------------------------
 */
void
Sw_HidDiscRV::tieParametersWith(RV*const other,
				bool checkStructure)
{
  assert ( other -> discrete() );
  assert ( other -> switching() );
  assert ( !other -> hidden() );

  if (!SwDiscRV::tieParametersWith( (SwDiscRV*)other,checkStructure))
    error("ERROR: trying to tie parameters of Switching Hidden Discrete RV '%s' with RV '%s' but they have different structure.",
	  name().c_str(),other->name().c_str());
  curCPT = ((Sw_HidDiscRV*)other)->curCPT;
}

#endif


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
Sw_HidDiscRV* Sw_HidDiscRV::cloneRVShell()
{
  Sw_HidDiscRV* rv = (Sw_HidDiscRV*)HidDiscRV::cloneRVShell();
  rv->tieParametersWith(this);
  return rv;
}

