/*-
 * GMTK_Sparse1DPMF.cc
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
 * Sparse1DPMF::Sparse1DPMF()
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
Sparse1DPMF::Sparse1DPMF() 
{

}



/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::read(is)
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
Sparse1DPMF::read(iDataStreamFile& is)
{

  is.read(length,"Sparse1DPMF::read, distribution length");

  if (length <= 0)
    error("Sparse1DPMF: read length (%d) < 0 in input",length);

  pmf.resize(length);

  for (int i=0;i<length;i++) {
    double val;
    is.readDouble(val,"Sparse1DPMF::read, reading value");
    if (val < 0 || val > 1)
      error("Sparse1DPMF: read, invalid pmf value (%g)",val);
    pmf[i] = val;
  }
}




/*-
 *-----------------------------------------------------------------------
 * Sparse1DPMF::read(is)
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
Sparse1DPMF::write(oDataStreamFile& os)
{
  assert (nFeats > 0);

  os.write(length,"Sparse1DPMF::write, distribution length");
  for (int i=0;i<length;i++) {
    // convert out of log domain and write out.
    os.writeDouble(pmf[i].unlog(),"Sparse1DPMF::write, writing value");
  }
  os.nl();

}





////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

