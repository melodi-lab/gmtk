/*-
 * GMTK_DoubleMatrix.cc
 *     General matrix class (for anything that needs generic scalars, vectors, or matrices).
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
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

#include "GMTK_DoubleMatrix.h"
#include "GMTK_GMParms.h"
#include "tieSupport.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif
VCID(HGID)



////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * DoubleMatrix::DoubleMatrix()
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
DoubleMatrix::DoubleMatrix() 
{
  numTimesShared = 0;
  refCount = 0;
}


/*-
 *-----------------------------------------------------------------------
 * DoubleMatrix::read(is)
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
DoubleMatrix::read(iDataStreamFile& is)
{
  NamedObject::read(is);
  is.read(_rows,"DoubleMatrix::read, distribution rows");
  if (_rows <= 0)
    error("DoubleMatrix: read rows (%d) < 0 in input",_rows);

  is.read(_cols,"DoubleMatrix::read, distribution cols");
  if (_cols <= 0)
    error("DoubleMatrix: read cols (%d) < 0 in input",_cols);

  values.resize(_rows*_cols);

  // use vector read
  is.read(values.ptr,values.len(),"DoubleMatrix::read, reading values");
  setBasicAllocatedBit();
  numTimesShared = 0;
  refCount = 0;
}



/*-
 *-----------------------------------------------------------------------
 * DoubleMatrix::write(os)
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
DoubleMatrix::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  NamedObject::write(os);
  os.write(_rows,"DoubleMatrix::write, distribution rows");
  os.write(_cols,"DoubleMatrix::write, distribution cols");

  // os.write(arr.ptr,arr.len(),"DoubleMatrix::write, writeing value");

  double * ptr  = values.ptr;
  for (int i=0; i < _rows; i++) {
    os.write(ptr,_cols,"DoubleMatrix: writing a row");
    ptr += _cols;
    os.nl();
  }
}



////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////


/*-
 *-----------------------------------------------------------------------
 * DoubleMatrix::cleanClone()
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
DoubleMatrix*
DoubleMatrix::cleanClone()
{
  assert ( basicAllocatedBitIsSet() );

  // TODO: when cloning is working, need to modify MixtureCommon to
  // keep track of when this object has already been cloned in a
  // training iteration.

  DoubleMatrix* clone = new DoubleMatrix();
  setName(new_name(name(),&GM_Parms.doubleMatsMap));
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
/// for this object. Since different object users of a double matrix might be very different
//  there is nothign more that is done here other than the most generic of EM support
//  (to ensure that the bits are set correctly) and we assume that the user will do
//  anything here further from the outside.


void DoubleMatrix::emStartIteration() {
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
DoubleMatrix::emIncrement(logpr prob) {
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emOnGoingBitIsSet())
    emStartIteration();
  accumulatedProbability += prob;
}

void DoubleMatrix::emEndIteration() 
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

void DoubleMatrix::emSwapCurAndNew() 
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  if (!emSwappableBitIsSet())
    return;
  values.swapPtrs(nextValues);
  emClearSwappableBit();
}
