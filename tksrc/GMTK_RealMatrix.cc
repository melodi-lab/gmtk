/*-
 * GMTK_RealMatrix.cc
 *     General matrix class (for anything that needs generic scalars, vectors, or matrices).
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
#include "GMTK_GMParms.h"
#include "tieSupport.h"

VCID("$Header$")


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
  numTimesShared = 0;
  refCount = 0;
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
  NamedObject::read(is);
  is.read(_rows,"RealMatrix::read, distribution rows");
  if (_rows <= 0)
    error("RealMatrix: read rows (%d) < 0 in input",_rows);

  is.read(_cols,"RealMatrix::read, distribution cols");
  if (_cols <= 0)
    error("RealMatrix: read cols (%d) < 0 in input",_cols);

  values.resize(_rows*_cols);

  // use vector read
  is.read(values.ptr,values.len(),"RealMatrix::read, reading values");
  setBasicAllocatedBit();
  numTimesShared = 0;
  refCount = 0;
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
  assert ( basicAllocatedBitIsSet() );

  NamedObject::write(os);
  os.write(_rows,"RealMatrix::write, distribution rows");
  os.write(_cols,"RealMatrix::write, distribution cols");

  // os.write(arr.ptr,arr.len(),"RealMatrix::write, writeing value");

  float * ptr  = values.ptr;
  for (int i=0; i < _rows; i++) {
    os.write(ptr,_cols,"RealMatrix: writing a row");
    ptr += _cols;
    os.nl();
  }
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * RealMatrix::cleanClone()
 *      make an exact clone of this object
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
RealMatrix*
RealMatrix::cleanClone()
{
  assert ( basicAllocatedBitIsSet() );

  // TODO: when cloning is working, need to modify MixtureCommon to
  // keep track of when this object has already been cloned in a
  // training iteration.

  RealMatrix* clone = new RealMatrix();
  setName(new_name(name(),&GM_Parms.realMatsMap));
  clone->_rows = _rows;
  clone->_cols = _rows;
  clone->refCount = 0;
  clone->numTimesShared = 0;
  clone->values.copyOtherIntoSelf(values);
  clone->setBasicAllocatedBit();

  // also add self to GMParms object.
  GM_Parms.add(clone);

  return clone;
}

/////////////////////////////////////////////////////////////////////////
/// EM Support: the next set of routines basically provide generic EM support
/// for this object. Since different object users of a real matrix might be very different
//  there is nothign more that is done here other than the most generic of EM support
//  (to ensure that the bits are set correctly) and we assume that the user will do
//  anything here further from the outside.


void RealMatrix::emStartIteration() {
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet()) return;
  if (emOnGoingBitIsSet()) return;
  if (!emEmAllocatedBitIsSet()) {
    nextValues.resize(values.size());
    emSetEmAllocatedBit();
  }
  emSetOnGoingBit();
  emSetSwappableBit();
  accumulatedProbability = 0.0;  
}

void
RealMatrix::emIncrement(logpr prob) {
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emOnGoingBitIsSet())
    emStartIteration();
  accumulatedProbability += prob;
}

void RealMatrix::emEndIteration() 
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emOnGoingBitIsSet())
    return;
  accumulatedProbability.floor();
  // stop EM
  emClearOnGoingBit();
}

void RealMatrix::emSwapCurAndNew() 
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emSwappableBitIsSet())
    return;
  values.swapPtrs(nextValues);
  emClearSwappableBit();
}
