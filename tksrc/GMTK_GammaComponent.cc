/*-
 * GMTK_GammaComponent.cc
 *        Any of the common code for the family of Gamma-like
 *        components classes.
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
#include <ctype.h>

#include <string>

#include "general.h"
VCID("$Header$")
#include "error.h"
#include "rand.h"

#include "GMTK_GammaComponent.h"
#include "GMTK_GMParms.h"
#include "GMTK_MixtureCommon.h"
#include "tieSupport.h"


////////////////////////////////////////////////////
// The default value of the minimum possible variance of any
// Gamma Distribution. This must be >= FLT_MIN for numeric stability,
// and it should be made availalbe as a command line parameter at some point.
double GammaComponent::_varianceFloor = GammaComponent::setVarianceFloor(1e-10);

double GammaComponent::setVarianceFloor(const double floor) 
{ 
  if (floor < FLT_MIN) 
    _varianceFloor = FLT_MIN; 
  else 
    _varianceFloor = floor; 
  return _varianceFloor;
}

// TODO: these values should not be exported to the command line
// until (and if) sharing is done.
// true if when a clone occurs, we use the same shape (i.e.,
// share the shape and only clone other things).
bool GammaComponent::cloneShareShape = false;
// true if when a clone occurs, we use the same scale
// (i.e., share the scale and clone other things)
bool GammaComponent::cloneShareScale = false;


/*-
 *-----------------------------------------------------------------------
 * recomputeDenominators()
 *      Recompute pre-computed values for probabilit evaluation.
 * 
 * Preconditions:
 *      Basic shape and scale parameter matrices should be allocated and current.
 *
 * Postconditions:
 *      denominator is recomputed.
 *
 * Side Effects:
 *      modifies variable 'denominator'
 *
 * Results:
 *      none
 *
 *-----------------------------------------------------------------------
 */

void
GammaComponent::recomputeDenominators()
{
  assert ( basicAllocatedBitIsSet() );
  denominators.resizeIfDifferent(_dim);
  for (unsigned i=0; i< _dim; i++) {
    // compute log_e((theta^k gamma(k)));
    // = k*log_e(theta) + log_e(gamma(k))
    const double k = shape->values.ptr[i];
    const double theta = scale->values.ptr[i];
    denominators.ptr[i] =  lgamma(k) +  k * log(theta);
  }
}


/*-
 *-----------------------------------------------------------------------
 * log_p()
 *      Computes the probability of this Gamma distribution vector (multiplying them together)
 * 
 * Preconditions:
 *      none
 *
 * Postconditions:
 *      nil
 *
 * Side Effects:
 *      nil, other than possible FPEs if the values are garbage
 *
 * Results:
 *      Returns the probability.
 *
 *-----------------------------------------------------------------------
 */

logpr
GammaComponent::log_p(const float *const x,
		    const Data32* const base,
		    const int stride)
{
  assert ( basicAllocatedBitIsSet() );

  double val = 0;

  // cache value
  const double localLower = lower;

  for (unsigned i = 0; i < _dim ; i++ ) {

    // enforce strict positivity.
    if (x[i] <= localLower) {
      // unfortnately there is no way here to give a better error
      // message.  I suppose if we threw an exception, it could be
      // caught where more information is availale regarding the
      // observation location and file name, etc.
      error("ERROR: obsevation value is %f (element %d), but a Gamma distribution (%s) has lower limit %f, all values must be strictly greater.",x[i],i,_name.c_str(),localLower);
    }

    // the gamma function is 
    //    x^{k-1} exp(-x/theta)/(theta^k gamma(k)) = (x/theta)^{k-1}exp(- (x/theta))/gamma(k)
    // theta is the (obviously trival) scale parameter, and where k is the shape parameter.
    // We note that k = 1 represents the exponential density for various scales theta. Note,
    // sometimes this density is represented with alpha = 1/theta. 

    const double x_val = x[i] - localLower; 
    const double k = shape->values.ptr[i];
    const double theta = scale->values.ptr[i];  
    // Compute the probability in the log domain.
    val += (k-1)*::log(x_val)  - x_val/theta - denominators.ptr[i];
  }

  return logpr(0,val);
}


// Slightly perturb either the shape or scale parameters.
void
GammaComponent::perturbShape()
{
  assert ( basicAllocatedBitIsSet() );
  for (unsigned i=0; i< _dim ; i++) {
    shape->values.ptr[i] +=  (shape->values.ptr[i]/8.0)*(rnd.drand48()-0.5); 
  }  
}
void
GammaComponent::perturbScale()
{
  assert ( basicAllocatedBitIsSet() );
  for (unsigned i=0;i<_dim;i++) {
    scale->values.ptr[i] +=  (scale->values.ptr[i]/8.0)*(rnd.drand48()-0.5); 
  }
}




void
GammaComponent::read(iDataStreamFile& is)
{
  // The index in the global scale 
  int scaleIndex; 
  // The index in the global shape
  int shapeIndex;

  // read name
  NamedObject::read(is);

  // read range (lower,infty)
  is.read(lower,"beta lower");

  // read scale vector
  string str;
  is.read(str);


  if (GM_Parms.realMatsMap.find(str) ==  GM_Parms.realMatsMap.end()) 
      error("Error: GammaComponent '%s' in file '%s' line %d specifies scale name '%s' that does not exist",
	    _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  scaleIndex = GM_Parms.realMatsMap[str];
  scale = GM_Parms.realMats[scaleIndex];

  if (scale->rows() != 1 || (unsigned)scale->cols() != _dim)
    error("Error: GammaComponent '%s' in file '%s' line %d specifices a scale real matrix '%s' with rows,cols = %d, %d but rows must be 1 and columns must be %d.\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  scale->name().c_str(),
	  scale->rows(),
	  scale->cols(),
	  _dim);


  // read shape vector
  is.read(str);
  if (GM_Parms.realMatsMap.find(str) == GM_Parms.realMatsMap.end())
    error("Error: GammaComponent '%s' in file '%s' line %d specifies shape name '%s' that does not exist",
	  _name.c_str(),is.fileName(),is.lineNo(),str.c_str());
  
  shapeIndex = GM_Parms.realMatsMap[str];
  shape = GM_Parms.realMats[shapeIndex];

  if (shape->rows() != 1 || (unsigned)shape->cols() != _dim)
    error("Error: GammaComponent '%s' in file '%s' line %d specifices a shape real matrix '%s' with rows,cols = %d, %d but rows must be 1 and columns must be %d.\n",
	  _name.c_str(),is.fileName(),is.lineNo(),
	  shape->name().c_str(),
	  shape->rows(),
	  shape->cols(),
	  _dim);

  setBasicAllocatedBit();
  recomputeDenominators();
}


void
GammaComponent::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );

  // write the type of self and the name
  os.write((int)Component::GammaComponent);
  NamedObject::write(os);
  os.nl();
  os.write(lower);  os.nl();

  // write object names
  os.write(scale->name());
  os.write(shape->name());
  os.nl();
}


/*-
 *-----------------------------------------------------------------------
 * GammaComponent::makeRandom
 *      calls makeRandom for the parameters
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
GammaComponent::makeRandom()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;
  for (unsigned i = 0; i < _dim ; i++ ) {
    scale->values.ptr[i] = rnd.drand48pe();
    shape->values.ptr[i] = rnd.drand48pe();
  }
}


/*-
 *-----------------------------------------------------------------------
 * makeUniform()
 *      calls makeUniform for the parameters.
 * 
 * Preconditions:
 *      Object must be allocated.
 *
 * Postconditions:
 *      Object has uniform values.
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
GammaComponent::makeUniform()
{
  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  for (unsigned i = 0; i < _dim; i++ ) {
    // not really uniform, but a surrogate. 
    scale->values.ptr[i] = 2.0;
    shape->values.ptr[i] = 1.0;
  }

}


/*-
 *-----------------------------------------------------------------------
 * noisyClone()
 *      Create a copy of self, but with cloned perturbed the mean/variance vectors
 * 
 * Preconditions:
 *      Obj must be read in, and basicAllocatedBitIsSet() must be true.
 *
 * Postconditions:
 *      none.
 *
 * Side Effects:
 *      'this' is not changed at all. Allocates new memory though.
 *
 * Results:
 *      returns the new noisy mean.
 *
 *-----------------------------------------------------------------------
 */
Component*
GammaComponent::noisyClone()
{
  assert ( basicAllocatedBitIsSet() );

  map<Component*,Component*>::iterator it = MixtureCommon::mcCloneMap.find(this);

  // first check if self is already cloned, and if so, return that.
  if (it == MixtureCommon::mcCloneMap.end()) {
    GammaComponent* clone;
    clone = new GammaComponent(_dim);

    clone->lower = lower;

    unsigned cloneNo=0; do {
      char buff[256];
      sprintf(buff,"%d",cloneNo);
      clone->_name = _name + string("_cl") + buff;
      cloneNo++;
    } while (GM_Parms.componentsMap.find(clone->_name) 
	     != GM_Parms.componentsMap.end());

    // sharing during training is not implemented (yet). Make sure the
    // user doesn't change code above and set these values.
    assert (!cloneShareShape && !cloneShareScale);

    if (cloneShareShape && cloneShareScale) {
      warning("WARNING: Gamma Component '%s' is cloning, and was asked to share both scale and shape. No sharing is occuring instead.",name().c_str());
      clone->scale = scale->cleanClone();
      clone->shape = shape->cleanClone();
      clone->perturbScale();
      clone->perturbShape();
    } else {
      if (cloneShareShape)
	clone->shape = shape;
      else {
	clone->shape = shape->cleanClone();
	clone->perturbShape();
      }

      if (cloneShareScale)
	clone->scale = scale;
      else {
	clone->scale = scale->cleanClone();
	clone->perturbScale();
      }
    }
    
    // need to tell mean, and covar that either
    //    1) if this is a new mean,covar that a
    //       a parent object is using them, or
    //    2) if this is a shared mean,covar that
    //       an additional parent object is using them.
    clone->scale->numTimesShared++;
    clone->shape->numTimesShared++;

    clone->setBasicAllocatedBit();
    MixtureCommon::mcCloneMap[this] = clone;

    // also add self to GMParms object.
    GM_Parms.add(clone);

    return clone;
  } else {
    return (*it).second;
  }
}



/*-
 *-----------------------------------------------------------------------
 * GammaComponent::identicalIndependentClone
 *      creates an exact copy of this object that shares nothing with
 *      the original
 *
 * Preconditions:
 *      1) object being copied should be allocated
 *      2) GM_Parms should contain all parameters, so that a unique name
 *         for the new object can be generated
 *
 * Postconditions:
 *      none
 *
 * Side Effects:
 *      the new object is added to GM_Parms
 *
 * Results:
 *      a pointer the new object
 *
 *-----------------------------------------------------------------------
 */
Component* 
GammaComponent::identicalIndependentClone()
{
  GammaComponent* newDG = new GammaComponent(_dim);

  newDG->shape = shape->cleanClone();
  newDG->scale = scale->cleanClone();
  newDG->lower = lower;

  // don't change usage counts here - do it in calling function,
  // because only that knows how the sharing is arranged

  newDG->_name = new_name(name(),&GM_Parms.componentsMap);
  newDG->setBasicAllocatedBit();
  GM_Parms.add(newDG);

  return newDG;
}

/////////////////
// EM routines //
/////////////////



void
GammaComponent::emStartIteration()
{

  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (emOnGoingBitIsSet())
    return;

  if (!emEmAllocatedBitIsSet()) {
    // this is presumably the first time
    emSetEmAllocatedBit();
  }

  // EM iteration is now going.
  emSetOnGoingBit();
  emSetSwappableBit();

  accumulatedProbability = 0.0;

  sumx.resizeIfDifferent(_dim);
  sumxx.resizeIfDifferent(_dim);
  sumlogx.resizeIfDifferent(_dim);
  for (unsigned i = 0 ; i < _dim; i ++ ) {
    sumx.ptr[i] = sumxx.ptr[i] = sumlogx.ptr[i] = 0;

  }

  scale->emStartIteration();
  shape->emStartIteration();

}


void
GammaComponent::emIncrement(logpr prob,
			  const float*f,
			  const Data32* const base,
			  const int stride)
{

  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  if (!emOnGoingBitIsSet())
    emStartIteration();

  if (prob < minIncrementProbabilty) {
    missedIncrementCount++;
    return;
    // don't accumulate anything since this one is so small and
    // if we did an unlog() and converted to a single precision
    // floating point number, it would be a denomral.
  }

  accumulatedProbability += prob;
  // prob.unlog() here so it doesn't need to be done
  // twice by the callees.
  const double fprob = prob.unlog();

  // we next do an estimate of the stats for a Gamma distribuion.
  // NOTE: this is currently just a hack, this is *NOT* the
  // MLE. To do an MLE of a Gamma requires an additional iterative
  // procedure.

  // cache value
  const double localLower = lower;

  if (!fisherKernelMode) {
    // do normal EM increment
    // this call is optimized to do the 1st and 2nd moment stats simultaneously
    for (unsigned i = 0; i < _dim; i++ ) {
      // enforce strict positivity.
      if (f[i] <= localLower) {
	// unfortnately there is no way here to give a better error
	// message.  I suppose if we threw an exception, it could be
	// caught where more information is availale regarding the
	// observation location and file name, etc.
	error("ERROR: em training, obsevation value is %f (element %d), but Gamma distribution (%s) has lower limit %f, all values must be strictly greater.",f[i],i,_name.c_str(),localLower);
      }

      const double x_val = f[i] - localLower; 
      sumx.ptr[i] += x_val*fprob;
      sumxx.ptr[i] += x_val*x_val*fprob;
      sumlogx.ptr[i] += ::log(x_val)*fprob;
    }
  } else {
    error("ERROR: fisher kernel not yet implemented with gamma components");
  }

  scale->emIncrement(prob);
  shape->emIncrement(prob);
}

void
GammaComponent::emEndIteration()
{

  assert ( basicAllocatedBitIsSet() );
  if (!emAmTrainingBitIsSet())
    return;

  // if EM not ongoing, we just return
  // here since, if this object is shared,
  // we don't want to end the epoch > 1 time
  if (!emOnGoingBitIsSet())
    return; 

  accumulatedProbability.floor();
  if (accumulatedProbability.zero()) {
    // Note: we assume here that the mean and covar object will check for us that 
    // its accumulated probability is above threshold. Here, we just
    // check for zero, and then issue a warning if needed. This might not
    // indicate a problem as the mean and covar of this object might be shared
    // and might have received plenty of count.
    warning("WARNING: Gamma Component named '%s' did not receive any accumulated probability in EM iteration. Global missed increment count is %d. Also check child scale '%s' and shape '%s'",
	    name().c_str(),
	    missedIncrementCount,
	    scale->name().c_str(),
	    shape->name().c_str());
  }

  // TODO: implement the sharing case properly.
  const double inv_denom = accumulatedProbability.inverse().unlog();

  for (unsigned i = 0; i < _dim; i++ ) {
    const double mean = sumx.ptr[i]*inv_denom;
    double variance = sumxx.ptr[i]*inv_denom - mean*mean;
    if (variance < _varianceFloor)
      variance = _varianceFloor;
    double meanlogx = sumlogx.ptr[i]*inv_denom;

    // A quick but hacky way to get the parameters is to do the
    // following:
    // 
    //    shape->nextValues.ptr[i] = mean*mean/variance;
    //    scale->nextValues.ptr[i] = variance/mean; // = mean/shape
    // 
    // A better way to get MLEs is to use a generalized Newton method
    // as suggested by Tom Minka.
    // 
    // We start with initial estimates based on the sample mean and
    // variance.
    double shapev = mean*mean/variance;
    // printf("shape = %f\n",shapev);
    double next_shapev;
    // do no more than this many iterations (in practice, it seems to do about 3 or 4 iters).
    unsigned max_iters = 10; 
    do {
      int dummy = 0;
      // do a Newton update.
      next_shapev = 1.0/(1/shapev + (meanlogx - log(mean) + log(shapev) - digamma(shapev,&dummy))/( shapev - shapev*shapev*trigamma(shapev,&dummy)));
      // printf("next_shape = %f\n",next_shapev);
      if (100.0*fabs(next_shapev - shapev)/shapev < 0.001)
	break;
      shapev = next_shapev;
    } while (max_iters--);

    // finalize he values.
    shape->nextValues.ptr[i] = shapev;
    scale->nextValues.ptr[i] = mean/shapev;

  }

  scale->emEndIteration();
  shape->emEndIteration();
  emClearOnGoingBit();

}


void
GammaComponent::emSwapCurAndNew()
{
  assert ( basicAllocatedBitIsSet() );

  if (!emAmTrainingBitIsSet())
    return;

  if (!emSwappableBitIsSet())
    return;

  scale->emSwapCurAndNew();
  shape->emSwapCurAndNew();
  recomputeDenominators();

  emClearSwappableBit();
}


////////////////////////////////////////////////////////////
// Parallel EM support
////////////////////////////////////////////////////////////


void
GammaComponent::emStoreObjectsAccumulators(oDataStreamFile& ofile,
					 bool writeLogVals,
					 bool writeZeros)
{
  assert (emEmAllocatedBitIsSet());
    

  // since this is a Gamma, we ignore the writeLogVals
  // argument since it doesn't make sense to take log of
  // these values since they are continuous. etc.
  if (writeZeros) {
    for (unsigned i=0;i<3*_dim;i++) {
      ofile.write(0.0,"Gamma Component store accum.");
    }
  } else {
    for (unsigned i = 0; i < _dim; i++) {
      ofile.write(sumx.ptr[i],"Gamma Component store accum.");
      ofile.write(sumxx.ptr[i],"Gamma Component store accum.");
      ofile.write(sumlogx.ptr[i],"Gamma Component store accum.");
    }
  }
}


void
GammaComponent::emLoadObjectsDummyAccumulators(iDataStreamFile& ifile)
{

  for (unsigned i = 0; i < _dim; i++) {
    double tmp;
    ifile.read(tmp,"Gamma load accums nm.");
    ifile.read(tmp,"Gamma load accums nc.");
    ifile.read(tmp,"Gamma load accums lg.");
  }
}


void
GammaComponent::emZeroOutObjectsAccumulators()
{
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    sumx.ptr[i] = sumxx.ptr[i] = sumlogx.ptr[i] = 0.0;
  }
}

void
GammaComponent::emLoadObjectsAccumulators(iDataStreamFile& ifile)
{
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    ifile.read(sumx.ptr[i],"Gamma Component load accums x");
    ifile.read(sumxx.ptr[i],"Gamma Component load accums xx.");
    ifile.read(sumlogx.ptr[i],"Gamma Component load accums lgx.");
  }
}


void
GammaComponent::emAccumulateObjectsAccumulators(iDataStreamFile& ifile)
{
  //
  assert (emEmAllocatedBitIsSet());
  for (unsigned i = 0; i < _dim; i++) {
    double tmp;
    ifile.read(tmp,"Gamma component accumulate accums sumx.");
    sumx.ptr[i] += tmp;
    ifile.read(tmp,"Gamma Gaussian accumulate accums sumxx.");
    sumxx.ptr[i] += tmp;
    ifile.read(tmp,"Gamma Gaussian accumulate accums sumlogx.");
    sumlogx.ptr[i] += tmp;
  }
}


////////////////////////////////////////////////////////////
// Sample generation
////////////////////////////////////////////////////////////



void GammaComponent::sampleGenerate(float *const sample,
				  const Data32* const base,
				  const int stride)
{
  error("GammaComponent sample generate not implemented");
}



////////////////////////////////////////////////////////////
// internal derivatives of log Gamma functions
////////////////////////////////////////////////////////////

//****************************************************************************80
double digamma ( double x, int *ifault )
//****************************************************************************80
//
//  Purpose:
//
//    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Jose Bernardo
//    FORTRAN90 version by John Burkardt
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double DIGAMA, the value of the digamma function at X.
//
{
  double c = 8.5;
  double d1 = -0.5772156649;
  double r;
  double s = 0.00001;
  double s3 = 0.08333333333;
  double s4 = 0.0083333333333;
  double s5 = 0.003968253968;
  double value;
  double y;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    *ifault = 1;
    return value;
  }
//
//  Initialize.
//
  *ifault = 0;
  y = x;
  value = 0.0;
//
//  Use approximation if argument <= S.
//
  if ( y <= s )
  {
    value = d1 - 1.0 / y;
    return value;
  }
//
//  Reduce to DIGAMA(X + N) where (X + N) >= C.
//
  while ( y < c )
  {
    value = value - 1.0 / y;
    y = y + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion if argument > C.
//
  r = 1.0 / y;
  value = value + log ( y ) - 0.5 * r;
  r = r * r;
  value = value - r * ( s3 - r * ( s4 - r * s5 ) );

  return value;
}
//****************************************************************************80


//****************************************************************************
double trigamma ( double x, int *ifault )
//****************************************************************************
//
//  Purpose:
//
//    TRIGAM calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
//
//  Modified:
//
//    19 January 2008
//
//  Author:
//
//    BE Schneider
//    Modifications by John Burkardt
//
//  Reference:
//
//    BE Schneider,
//    Algorithm AS 121:
//    Trigamma Function,
//    Applied Statistics, 
//    Volume 27, Number 1, pages 97-99, 1978.
//
//  Parameters:
//
//    Input, double X, the argument of the trigamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double TRIGAM, the value of the trigamma function at X.
//
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381;
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  z = x;
//
//  Use small value approximation if X <= A.
//
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
//
//  Increase argument to ( X + I ) >= B.
//
  value = 0.0;

  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
//
//  Apply asymptotic formula if argument is B or greater.
//
  y = 1.0 / z / z;

  value = value + 0.5 * 
      y + ( 1.0 
    + y * ( b2  
    + y * ( b4  
    + y * ( b6  
    + y *   b8 )))) / z;

  return value;
}
//****************************************************************************80
