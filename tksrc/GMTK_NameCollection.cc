/*-
 * GMTK_NameCollection.cc
 *     named colleciton.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
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



#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_NameCollection.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "GMTK_Mixture.h"


VCID("$Header$");

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * NameCollection::NameCollection()
 *      Constructor
 *
 * Results:
 *      Constructs the object.
 *
 * Side Effects:
 *      None so far.
 *
 *-----------------------------------------------------------------------
 */
NameCollection::NameCollection()
{
}


////////////////////////////////////////////////////////////////////
//        I/O
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * NameCollection::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the internal table
 *
 *-----------------------------------------------------------------------
 */

void
NameCollection::read(iDataStreamFile& is)
{

  NamedObject::read(is);
  int length;
  is.read(length,"NameCollection::read length");

  if (length <= 0) 
    error("ERROR: reading file '%s', NameCollection '%s', must have positive number of entries.",
	  is.fileName(),name().c_str(),length);
  table.resize(length);
  for (int i=0;i<length;i++) {
    is.read(table[i],"MTCPT::read table entry");
  }
}


/*-
 *-----------------------------------------------------------------------
 * NameCollection::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(table.size(),"NameCollection::write length");
  for (unsigned i=0;i<table.size();i++) {
    os.write(table[i],"NameCollection::write table entry");
    if ((i+1) % 10 == 0)
      os.nl();
  }
  os.nl();
}




/*-
 *-----------------------------------------------------------------------
 * NameCollection::fillMxTable
 *      fill the spmf table with entries from GM_Params
 *      This routine must be called to fill the table
 *      that allow a RV to go from an integer index value
 *      to a pointer to the appropriate object.
 *      Note that if mxTable is already filled (size greater 
 *      or equal to table), then nothing happens.
 *
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      fills the mxTable entry.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::fillMxTable()
{
  if (mxTable.size() >= table.size())
    return;

  mxTable.resize(table.size());
  for (unsigned i=0;i<table.size();i++ ) {
    GMParms::ObjectMapType::iterator it;
    if ((it = GM_Parms.mixturesMap.find(table[i])) == GM_Parms.mixturesMap.end())
      error("Error: collection '%s' has named a mixture '%s' (at table entry %d) that doesn't exist.",
	    name().c_str(),table[i].c_str(),i);
    unsigned indx = (*it).second;
    mxTable[i] = GM_Parms.mixtures[indx];
  }
}



/*-
 *-----------------------------------------------------------------------
 * NameCollection::fillspmfTable
 *      fill the mx table with entries from GM_Params.
 *      This routine must be called to fill the table
 *      that allow a RV to go from an integer index value
 *      to a pointer to the appropriate object.
 *      Note that if mxTable is already filled (size greater 
 *      or equal to table), then nothing happens.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      fills the mxTable entry.
 *
 *-----------------------------------------------------------------------
 */
void
NameCollection::fillSpmfTable()
{
  if (spmfTable.size() >= table.size())
    return;

  spmfTable.resize(table.size());
  for (unsigned i=0;i<table.size();i++ ) {
    GMParms::ObjectMapType::iterator it;
    if ((it = GM_Parms.sPmfsMap.find(table[i])) == GM_Parms.sPmfsMap.end())
      error("Error: collection '%s' has named an SPMF '%s' (at table entry %d) that doesn't exist.",
	    name().c_str(),table[i].c_str(),i);
    unsigned indx = (*it).second;
    spmfTable[i] = GM_Parms.sPmfs[indx];
  }
}

