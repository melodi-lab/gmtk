/*-
 * GMTK_Discrete1DPDF.cc
 *     Trainable (with say EM) 1D discrete probability
 *     distributions.
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


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Discrete1DPDF::Discrete1DPDF()
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
Discrete1DPDF::Discrete1DPDF() 
{

}



/*-
 *-----------------------------------------------------------------------
 * Discrete1DPDF::read(is)
 *      read in a distribution from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted into the log domain (so discrete data on disk
 *      is NOT represented as log probabilities. This is because 
 *          1) discrete data typically doesn't need such a huge dynamic range
 *              (like Gaussian probabilties do).
 *          2) it is easier to examine data on disk when it is not in log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the pmf member function in the object.
 *
 *-----------------------------------------------------------------------
 */
void
Discrete1DPDF::read(iDataStreamFile& is)
{

  is.read(length,"Discrete1DPDF::read, distribution length");

  if (length <= 0)
    error("Discrete1DPDF: read length (%d) < 0 in input",length);

  pmf.resize(length);

  for (int i=0;i<length;i++) {
    double val;
    is.readDouble(val,"Discrete1DPDF::read, reading value");
    if (val < 0 || val > 1)
      error("Discrete1DPDF: read, invalid pmf value (%g)",val);
    pmf[i] = val;
  }
}




/*-
 *-----------------------------------------------------------------------
 * Discrete1DPDF::read(is)
 *      write out distribution to file 'os'. 
 *      the data probs are stored on disk as doubles,  NOT in log domain.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effectcs other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
Discrete1DPDF::write(oDataStreamFile& os)
{
  assert (nFeats > 0);

  os.write(length,"Discrete1DPDF::write, distribution length");
  for (int i=0;i<length;i++) {
    // convert out of log domain and write out.
    os.writeDouble(pmf[i].unlog(),"Discrete1DPDF::write, writing value");
  }
  os.nl();

}





////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

