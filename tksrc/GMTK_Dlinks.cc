/*-
 * GMTK_DLINKS.cc
 *        2D arrays of 2-tuples containing <time lag, feature offset>
 *        The time lag says where, relative to the current position
 *        the feature dependency is, and the feature offset
 *        says to which feature relative to feature 0.
 *
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

#include "GMTK_Dlinks.h"
#include "GMTK_DlinkMatrix.h"

VCID("$Header$");


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * Dlinks::Dlinks()
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
Dlinks::Dlinks()
{

}



/*-
 *-----------------------------------------------------------------------
 * Dlinks::read(is)
 *      read in a distribution from file 'is'. 
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
Dlinks::read(iDataStreamFile& is)
{
  NamedObject::read(is);
  int nFeats;
  is.read(nFeats,"Dlinks::read, num feats");
  if (nFeats <= 0)
    error("Dlinks::read, read num feats (%d) < 0 in input",nFeats);

  dIndices.resize(nFeats);
  
  _minLag = 1000000;
  _maxLag = -1000000;
  for (int i=0;i<nFeats;i++) {
    int nLinks;
    is.read(nLinks,"Dlinks::read, nLinks");

    // Note we explicitely allow for there to be 0 links here.
    // If so, the array size will be set to have zero length.
    if (nLinks < 0)
      error("Dlinks::read, read nLinks (%d) < 0 in input",nLinks);
    dIndices[i].resize(nLinks);

    for (int j=0;j<nLinks;j++) {
      int l,o;
      // lags can be pos or negative
      is.read(l,"Dlinks::read, lag");      
      if (l < _minLag)
	_minLag = l;
      if (l > _maxLag)
	_maxLag = l;

      // offsets must be >= 0
      is.read(o,"Dlinks::read, offset");
      if (o < 0)
	error("Dlinks::read, read offset (%d) < 0 in input",o);
      dIndices[i][j].lag = l;
      dIndices[i][j].offset = o;
    }
  }
}



/*-
 *-----------------------------------------------------------------------
 * Dlinks::write(os)
 *      write out distribution to file 'os'. 
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
Dlinks::write(oDataStreamFile& os)
{
  NamedObject::write(os);
  os.write(dim(),"Dlinks::write, num feats");
  os.nl();
  for (int i=0;i<dim();i++) {
    os.write(numLinks(i),"Dlinks::write, nLinks");
    for (int j=0;j<numLinks(i);j++) {
      os.write(dIndices[i][j].lag,"Dlinks::write, lag");      
      os.write(dIndices[i][j].offset,"Dlinks::write, offset");
    }
    os.nl();
  }
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
Dlinks::compatibleWith(DlinkMatrix& d)
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
 * preCompute()
 *      precompute the offset array.
 * 
 * Preconditions:
 *      basic object should be read in.
 *
 * Postconditions:
 *      offset array is computed.
 *
 * Side Effects:
 *      destroys old offset array
 *
 * Results:
 *      retunrs nil
 *
 *-----------------------------------------------------------------------
 */
void
Dlinks::preCompute(const unsigned stride)
{
  // first go through and find out how long it needs to be
  unsigned len = 0;
  for (int i=0;i<dim();i++) 
    len += numLinks(i);

  preComputedOffsets.resize(len);

  unsigned entry = 0;
  for (int i=0;i<dim();i++)
    for (int j=0;j<numLinks(i);j++) {
      preComputedOffsets[entry++] = 
	dIndices[i][j].lag*stride+dIndices[i][j].offset;
    }
  assert ( len == entry );

}
