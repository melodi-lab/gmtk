/*-
 * GMTK_RealArray.cc
 *     General real array (for means, etc.)
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
#include <float.h>
#include <assert.h>

#include "general.h"
#include "error.h"

#include "GMTK_RealArray.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * RealArray::RealArray()
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
RealArray::RealArray() 
{


}


/*-
 *-----------------------------------------------------------------------
 * RealArray::read(is)
 *      read in the array from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted to the log domain.
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
RealArray::read(iDataStreamFile& is)
{
  int length;

  is.read(length,"RealArray::read, distribution length");
  if (length <= 0)
    error("RealArray: read length (%d) < 0 in input",length);

  resize(length);

  for (int i=0;i<length;i++) {
    is.read(operator[](i),"RealArray::read, reading value");
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

  os.write(len(),"RealArray::write, length");
  for (int i=0;i<len();i++) {
    os.write(operator[](i),"RealArray::write, values");
  }
  os.nl();
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
