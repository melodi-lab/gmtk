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

  if (strIsInt(str.c_str(),meanIndex)) {
    if (meanIndex < 0 || meanIndex >= (int)GM_Parms.means.size()) {
      error("Error: DiagGaussian '%s' specifies mean index (%d) that does not exist",
	    _name.c_str(),meanIndex);
    }
  } else {
    GMParms::MeansMapType::iterator it;
    it = GM_Parms.meansMap.find(str);
    if (it == GM_Parms.meansMap.end()) {
      error("Error: DiagGaussian '%s' specifies mean name '%s' that does not exist",
	    _name.c_str(),str.c_str());
    }
    meanIndex = (*it).second;
  }
  means = GM_Parms.means[meanIndex];


  // read covariance vector
  is.read(str);

  if (strIsInt(str.c_str(),covarIndex)) {
    if (covarIndex < 0 || covarIndex >= (int)GM_Parms.covars.size()) {
      error("Error: DiagGaussian '%s' specifies covar index (%d) that does not exist",
	    _name.c_str(),covarIndex);
    }
  } else {
    GMParms::CovarsMapType::iterator it;
    it = GM_Parms.covarsMap.find(str);
    if (it == GM_Parms.covarsMap.end()) {
      error("Error: DiagGaussian '%s' specifies covar name '%s' that does not exist",
	    _name.c_str(),str.c_str());
    }
    covarIndex = (*it).second;
  }
  covar = GM_Parms.covars[covarIndex];
}


void
DiagGaussian::write(oDataStreamFile& os)
{
  error("DiagGaussian::write not implemented");
}


void
DiagGaussian::preCompute()
{
}



void
DiagGaussian::makeRandom()
{
}


void
DiagGaussian::makeUniform()
{
}





//
// compute the log probability of x with stride 'stride'
// 
logpr
DiagGaussian::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  return 0.0;
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









