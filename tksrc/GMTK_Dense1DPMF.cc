/*-
 * GMTK_1D_Dist.cc
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




void
Discrete1DPDF::read(iDataStreamFile& is)
{
  assert (nFeats > 0);


  is.read(length);

  if (length <= 0) {
    error("Discrete1DPDF: read a lenght (%d) of < 0 in file (%s)",
	  length,is.fileName);
  }


  numComs = new int[nFeats];
  bIndices = new bindex*[nFeats];

  is.readInt(origNumComponents,"GaussianMixture::read ncmps");
  if (origNumComponents <= 0)
    error("GaussianMixture::read, must have at least 1 component in a mixture.");
  newNumComponents = curNumComponents = origNumComponents;
  mixCoeffVanishThreshold = 
    logpr((double)1.0/curNumComponents) /
    logpr(mixCoeffVanishRatio);

  mixCoeffs = new logpr[origNumComponents];
  components.resize(origNumComponents);
  for (int i=0;i<origNumComponents;i++)
    components[i] = new DiagGaussian();

  int i;

  _minLag = MAXINT;
  _maxLag = -MAXINT;
  _maxOffset = -MAXINT;
  sum_nComsp1 = 0;
  sum_nComsp1SqH = 0;
  for (i=0;i<nFeats;i++) {
    is.readInt(numComs[i],"GaussianMixture::read ncms");

    sum_nComsp1 += numComs[i]+1;
    sum_nComsp1SqH += (numComs[i]+1)*(numComs[i]+2)/2;

    bIndices[i] = new bindex[numComs[i]];
    for (int j=0;j<numComs[i];j++) {
      int l,o;
      is.readInt(l,"GaussianMixture::read lag");
      is.readInt(o,"GaussianMixture::read offset");
      
      if (l < _minLag)
	_minLag = l;
      if (l > _maxLag)
	_maxLag = l;
      if (o > _maxOffset)
	_maxOffset = o;

      bIndices[i][j].lag = l;
      bIndices[i][j].offset = o;
    }
  }
  if (_minLag == MAXINT)
    _minLag = 0;
  if (_maxLag == -MAXINT)
    _maxLag = 0;
  if (_maxOffset == -MAXINT)
    _maxOffset = 0;
  
  for (i=0;i<origNumComponents;i++) {
    int dummy;
    is.readInt(dummy,"GaussianMixture::read cmp #");
    if (dummy != i) {
      error("Error in BMM file, component number, got %d, expected %d",dummy,i);
    }
    double mxcoef;
    is.readDouble(mxcoef,"GaussianMixture::read alpha");
    mixCoeffs[i] = mxcoef;
    components[i]->setParent(this);
    components[i]->read(is);
  }

  z_cache.resize(sum_nComsp1 - nFeats);

  bitmask |= bm_basicAllocated;
}






////////////////////////////////////////////////////////////////////
//        Misc Support
////////////////////////////////////////////////////////////////////

