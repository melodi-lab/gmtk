/*-
 * GMTK_MixGaussians.cc
 *        Code for mixtures of Gaussians of a variety of types.
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
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"

#include "GMTK_MixGaussians.h"
#include "GMTK_GMParms.h"



void
MixGaussians::read(iDataStreamFile& is)
{
  // read name
  NamedObject::read(is);
}

void
MixGaussians::write(oDataStreamFile& os)
{
  error("MixGaussians::write not implemented");
  NamedObject::write(os);
}


void
MixGaussians::preCompute()
{
}


//
// compute the log probability of x with stride 'stride'
// 
logpr
MixGaussians::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  return 0.0;
}


/////////////////
// EM routines //
/////////////////


void
MixGaussians::emStartIteration()
{
  error("not implemented");
}


void
MixGaussians::emSwapCurAndNew()
{
  error("not implemented");
}


void
MixGaussians::emIncrement(RandomVariable* rv,
			  logpr prob)
{
  error("not implemented");
}


void
MixGaussians::emEndIteration()
{
  error("not implemented");
}


void
MixGaussians::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
MixGaussians::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
MixGaussians::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void MixGaussians::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}









