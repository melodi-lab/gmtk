/*-
 * GMTK_DlinkMatrix.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, make no representations about
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
#include "rand.h"

#include "GMTK_DlinkMatrix.h"
#include "GMTK_Dlinks.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Dlinkmatrix::Dlinkmatrix()
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
DlinkMatrix::DlinkMatrix() 
{
}




////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * compatibleWith()
 *      returns true of this object is compatible with the argument function.
 * 
 * Preconditions:
 *      Object must be read in.
 *
 * Postconditions:
 *      --
 *
 * Side Effects:
 *      none
 *
 * Results:
 *      true only if compatibility holds.
 *
 *-----------------------------------------------------------------------
 */
bool 
DlinkMatrix::compatibleWith(Dlinks& d)
{
  if (numFeats() != d.numFeats())
    return false;
  for (int i=0;i<numFeats();i++) {
    if (numLinks(i) != d.numLinks(i))
      return false;
  }
  return true;
}


/*-
 *-----------------------------------------------------------------------
 * makeRandom()
 *      assign random values to all elements
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has random values.
 *
 * Side Effects:
 *      destroys previous values.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::makeRandom()
{
  for (int i=0;i<numFeats();i++)
    for (int j=0;j<numLinks(i);j++)
      mat.arr[i][j] = rnd.drand48pe();
}



/*-
 *-----------------------------------------------------------------------
 * makeUniform()
 *      assign uniform (i.e., in this case 0) values to all elements
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has random values.
 *
 * Side Effects:
 *      destroys previous values.
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::makeUniform()
{
  for (int i=0;i<numFeats();i++)
    for (int j=0;j<numLinks(i);j++)
      mat.arr[i][j] = 0.0;
}

