/*-
 * GMTK_MeanVector.cc
 *     A vector used for means of Gaussians.
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
 * MeanVector::MeanVector()
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
MeanVector::MeanVector()
{


}


/*-
 *-----------------------------------------------------------------------
 * MeanVector::read(is)
 *      read in the mean
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes the 'means' member function in the object.
 *
 *-----------------------------------------------------------------------
 */

//  void
//  MeanVector::read(iDataStreamFile& is)
//  {}


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
  assert (nFeats > 0);

  os.write(length,"RealArray::write, length");
  for (int i=0;i<length;i++) {
    os.write(arr[i],"RealArray::write, values");
  }
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////







