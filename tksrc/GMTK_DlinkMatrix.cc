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
 * read(is)
 *      read in the array from file 'is'. 
 *      The data probs are stored as doubles, but when they are read in
 *      they are converted to the log domain. Also, they are read
 *      in as a single array for speed reasons.
 *
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      object is read in.
 *
 * Side Effects:
 *      Changes the pmf member function in the object.
 *
 * Results:
 *      nil
 *
 *-----------------------------------------------------------------------
 */
void
DlinkMatrix::read(iDataStreamFile& is)
{
  NamedObject::read(is);

  int _dim;

  is.read(_dim,"DlinkMatrix::read, _dim");
  if (_dim <= 0)
    error("ERROR: reading DlinkMatrix '%s' from file '%s', dim (%d) must be positive",name().c_str(),is.fileName(),_dim);

  _numLinks.resize(_dim);

  for (int i=0;i<_dim;i++) {
    int nlinks;
    is.read(nlinks,"DlinkMatrix::read, nlinks");

    if (nlinks < 0) 
      error("ERROR: reading DlinkMatrix '%s' from file '%s', # dlinks (%d) must be >= 0",name().c_str(),is.fileName(),nlinks);

    _numLinks[i] = nlinks;

    int oldLen = arr.len();
    arr.resizeAndCopy(oldLen+nlinks);

    for (int j=0;j<nlinks;j++) {
      is.read(arr[oldLen+j],"DlinkMatrix::read, v");
    }
  }
  setBasicAllocatedBit();
}




/*-
 *-----------------------------------------------------------------------
 * DlinkMatrix::write(os)
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
DlinkMatrix::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(dim(),"DlinkMatrix::write, dim()");
  os.nl();
  int ptr = 0;
  for (int i=0;i<dim();i++) {
    os.write(_numLinks[i],"DlinkMatrix::write, nlinks");
    for (int j=0;j<_numLinks[i];j++) {
      os.write(arr[ptr++],"DlinkMatrix::write val ");
    }
    os.nl();
  }
}




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
  if (dim() != d.dim())
    return false;
  for (int i=0;i<dim();i++) {
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
  for (int i=0;i<arr.len();i++)
    arr[i] = rnd.drand48pe();
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
  for (int i=0;i<arr.len();i++)
    arr[i] = 0;
}

