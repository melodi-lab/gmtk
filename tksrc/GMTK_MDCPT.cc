/*-
 * GMTK_MDCPT.cc
 *     A Multi-Dimensional dense Conditional Probability Table class.
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



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ieeefp.h>
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"


VCID("$Header$");


int MDCPT::warningDimensionality = 50;

////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * MDCPT::MDCPT()
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
MDCPT::MDCPT()
{


}


/*-
 *-----------------------------------------------------------------------
 * MDCPT::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the 'mdcpt' member function in the object.
 *
 *-----------------------------------------------------------------------
 */

void
MDCPT::read(iDataStreamFile& is)
{

  is.read(dimensionality,"MDCPT::read dimensionality");

  if (dimensionality <= 0) 
    error("MDCPT: read, trying to use 0 or negative (%d) dimensional table.",dimensionality);
  if (dimensionality >= warningDimensionality)
    warning("MDCPT: read, creating a %d-dimensional array.",dimensionality);

  cardinalities.resize(dimensionality);

  // read the cardinalities
  int numValues = 1;
  for (int i=0;i<dimensionality;i++) {
    is.read(cardinalities[i],"MDCPT::read cardinality");
    if (cardinality[i] <= 0)
      error("MDCPT: read, trying to use 0 or negative (%d) cardinality table.",cardinality[i]);
    numValues *= cardinality[i];
  }

  // Finally read in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1 here, if
  // we didn't use a large 1D loop. 
  mdcpt.resize(numValues);
  for (int i=0;i<numValues;i++) {
    double val;
    is.readDouble(val,"MDCPT::read, reading value");
    if (val < 0 || val > 1)
      error("MDCPT: read, invalid pmf value (%g)",val);
    mdcpt[i] = val;
  }

}


/*-
 *-----------------------------------------------------------------------
 * RealArray::write(os)
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
RealArray::write(oDataStreamFile& os)
{
  os.write(dimensionality,"MDCPT::write dimensionality");
  for (int i=0;i<dimensionality;i++) {
    os.write(cardinalities[i],"MDCPT::write cardinality");
  }

  // Finally write in the probability values (stored as doubles).
  // NOTE: We could check that things sum to approximately 1 here, if
  // we didn't use a large 1D loop. 
  for (int i=0;i<mdcpt.len();i++) {
    os.writeDouble(mdcpt[i].unlog(),"MDCPT::write, writing value");
  }

}


////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////







