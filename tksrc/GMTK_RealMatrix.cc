/*-
 * GMTK_RealMatrix.cc
 *     General matrix class (for means, etc.)
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

#include "GMTK_RealMatrix.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * RealMatrix::RealMatrix()
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
RealMatrix::RealMatrix() 
{


}


/*-
 *-----------------------------------------------------------------------
 * RealMatrix::read(is)
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
RealMatrix::read(iDataStreamFile& is)
{
  is.read(rows,"RealMatrix::read, distribution rows");
  if (rows <= 0)
    error("RealMatrix: read rows (%d) < 0 in input",rows);

  is.read(cols,"RealMatrix::read, distribution cols");
  if (cols <= 0)
    error("RealMatrix: read cols (%d) < 0 in input",cols);

  arr.resize(rows*cols);

  float *ptr = arr.ptr;
  for (int i=0;i<rows*cols;i++) {
    float val;
    is.read(val,"RealMatrix::read, reading value");
    *ptr++ = val;
  }
}




/*-
 *-----------------------------------------------------------------------
 * RealMatrix::write(os)
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
RealMatrix::write(oDataStreamFile& os)
{

  os.write(rows,"RealMatrix::write, distribution rows");
  os.write(cols,"RealMatrix::write, distribution cols");

  for (int i=0;i<rows*cols;i++) {
    os.write(arr[i],"RealMatrix::write, writeing value");
  }
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////
