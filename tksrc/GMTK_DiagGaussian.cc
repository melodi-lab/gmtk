/*-
 * GMTK_DiagGaussian.cc
 *        Code for plain vanilla diagonal Gaussians.
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
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$");
#include "error.h"
#include "rand.h"

#include "GMTK_DiagGaussian.h"
#include "GMTK_GMParms.h"



void
DiagGaussian::read(iDataStreamFile& is)
{
  // read name
  NamedObject::read(is);

  // read mean vector
  string str;
  is.read(str);

  if (GM_Parms.meansMap.find(str) ==  GM_Parms.meansMap.end()) 
      error("Error: DiagGaussian '%s' specifies mean name '%s' that does not exist",
	    _name.c_str(),str.c_str());
  meanIndex = GM_Parms.meansMap[str];
  means = GM_Parms.means[meanIndex];


  // read covariance vector
  is.read(str);
  if (GM_Parms.covarsMap.find(str) == GM_Parms.covarsMap.end())
    error("Error: DiagGaussian '%s' specifies covar name '%s' that does not exist",
	  _name.c_str(),str.c_str());
  
  covarIndex = GM_Parms.covarsMap[str];
  covar = GM_Parms.covars[covarIndex];

  // check that lengths match, etc.
  if (covar->dim() != means->dim()) {
    error("Error: DiagGaussian '%s' specifices a mean '%s' with dim %d and covariance '%s' with dim '%d'\n",
	  _name.c_str(),
	  means->name().c_str(),
	  means->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }
  if (covar->dim() != _dim) {
    error("Error: DiagGaussian '%s' of dim %d does not match its mean '%s' with dim %d or covariance '%s' with dim '%d'\n",
	  _name.c_str(),
	  _dim,
	  means->name().c_str(),
	  means->dim(),
	  covar->name().c_str(),
	  covar->dim());
  }

}


void
DiagGaussian::write(oDataStreamFile& os)
{
  error("DiagGaussian::write not implemented");
}


void
DiagGaussian::preCompute()
{
  covar->preCompute();
}



void
DiagGaussian::makeRandom()
{
  means->makeRandom();
  covar->makeRandom();
}


void
DiagGaussian::makeUniform()
{
  means->makeUniform();
  covar->makeUniform();
}



// This can be changed from 'float' to 'double' to
// provide extra range for temporary accumulators. Alternatively,
// decreasing the program's mixCoeffVanishRatio at the beginning
// of training should eliminate any component that produces
// such low scores.
#define DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE double


//
// compute the log probability of x with stride 'stride'
// 
logpr
DiagGaussian::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  logpr rc;
  rc.set_to_zero();

  // covariances must be precomputed for this
  // to work.
  const float *xp = x;
  const float *const x_endp = x + _dim;
  const float *mean_p = means->basePtr();
  const float *var_inv_p = covar->baseVarInvPtr();
  DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE d=0.0;
  do {
    const DIAG_GAUSSIAN_TMP_ACCUMULATOR_TYPE tmp
      = (*xp - *mean_p);
    d += tmp*tmp*(*var_inv_p);

    xp++;
    mean_p++;
    var_inv_p++;
  } while (x != x_endp);
  d *= -0.5;
  return logpr(0,(covar->log_inv_normConst() + d));
}



/////////////////
// EM routines //
/////////////////


void
DiagGaussian::emStartIteration()
{
  error("not implemented");
}


void
DiagGaussian::emSwapCurAndNew()
{
  error("not implemented");
}


void
DiagGaussian::emIncrement(RandomVariable* rv,
			  logpr prob)
{
  error("not implemented");
}


void
DiagGaussian::emEndIteration()
{
  error("not implemented");
}


void
DiagGaussian::emStoreAccumulators(oDataStreamFile& ofile)
{
  error("not implemented");
}

void
DiagGaussian::emLoadAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


void
DiagGaussian::emAccumulateAccumulators(iDataStreamFile& ifile)
{
  error("not implemented");
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void DiagGaussian::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("not implemented");
}









