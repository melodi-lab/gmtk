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

  is.read(numParents,"MDCPT::read numParents");

  if (numParents < 0) 
    error("MDCPT: read, trying to use negative (%d) num parents.",numParents);
  if (numParents >= warningNumParents)
    warning("MDCPT: read, creating MDCPT with %d parents",numParents);

  cardinalities.resize(numParents+1);
  cumulativeCardinalities.resize(numParents);

  // read the cardinalities
  int numValues = 1;
  for (i=0;i<=numParents;i++) {
    is.read(cardinalities[i],"MDCPT::read cardinality");
    if (cardinality[i] <= 0)
      error("MDCPT: read, trying to use 0 or negative (%d) cardinality table.",cardinality[i]);
    numValues *= cardinality[i];
  }
  for (int i=(numParents-1);i>=0;i--)
    cumulativeCardinalities[i] = 
      cumulativeCardinalities[i+1]*cardinality[i]

      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
 * MDCPT::write(os)
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
MDCPT::write(oDataStreamFile& os)
{
  os.write(numParents,"MDCPT::write numParents");
  for (int i=0;i<=numParents;i++) {
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
//        Probability Evaluation
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * MDCPT::setParentValues()
 *      Adjusts the current structure so that subsequent calls of
 *      probability routines will be conditioned on the given
 *      assigment to parent values.
 *  
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the mdcpt_ptr
 *
 *-----------------------------------------------------------------------
 */
void
MDCPT::setParentValues( sArray<int>& parentValues)
{

  assert ( parentValues.len() == numParents );
  
  





}

////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////







