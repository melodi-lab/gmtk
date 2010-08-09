/*-
 * GMTK_GammaComponent
 *        Gamma Component elements that are shared by all Gamma type component clases.
 *
 *  Written by Jeff Bilmes <bilmes@ee.washington.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2008, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

/*
 * A note on sharing: If no training is going on, it is possible to set up any form
 * of sharing that is specifyable in the parameter files. I.e., since a Gamma distribution
 * involves two parameters, we could have two Gamma distributions that share a shape but have
 * unique scale parameters. Since no training is going on, however, the utility of doing this
 * is questionable.
 *
 * During training, however, some forms of sharing currently won't work. What will work is
 * if two separate random variables share the same Gamma distribution (this can be at different
 * times or positions in a structure file). What won't work, however, is if two Gamma distributions
 * share an shape and/or a scale parameter vector. In other words, each real matrix (which
 * are used for the parameters of a given Gamma distribution, should be used only one time
 * during parameter training).
 *
 * A TODO is to derive the update equations for the constrained optimization case of
 * a Gamma distrubtion where one (or more) of the parameters are shared. See
 * also the Beta observation distribution.
 */



#ifndef GMTK_GAMMACOMPONENT_H
#define GMTK_GAMMACOMPONENT_H

#include "fileParser.h"
#include "logp.h"

#include "machine-dependent.h"

#include "GMTK_NamedObject.h"
#include "GMTK_Component.h"
#include "GMTK_RealMatrix.h"

// Derivatives of digamma generalized are polygamma functions such as
// trigamma.  We need them here.
extern double digamma ( double x, int *ifault );
extern double trigamma ( double x, int *ifault );

class GammaComponent :  public Component {

  ///////////////////////////////////////////////////////
  // The value that, if any variances go below, cause
  // either warnings (and adjustments) or errors to occur
  // depending on how probable this component is.
  static double _varianceFloor;

  // TODO: shared parmaeters during training is not
  // yet implemented!!!!

  // The scale parameter, which might be tied with
  // other Gamma distributions. We use a general 1x1 sized
  // real matrix for this.
  RealMatrix* scale;

  // The shape parameter, which might be tied with
  // other Gamma distributions. We use a general 1x1 sized
  // weight matrix for this.
  RealMatrix* shape;

  // set the left-most open interval limit to this value.
  // I.e., rather than distribution from (0,infty) we have
  // a distribution from (lower,infty)
  double lower;

  // used for EM accumulation
  sArray <double> sumx;
  sArray <double> sumxx;  
  sArray <double> sumlogx;


  // precomputed log denominator for the Gamma distribution.
  sArray <double> denominators;
  void recomputeDenominators();

  // Slightly perturb either the shape or scale parameters.
  void perturbShape();
  void perturbScale();


  /////////////////////////////////////////////////
  // modify the usage counts of any members that use them; typically
  // called with amount=1 or -1
  void adjustNumTimesShared(int amount){
    scale->numTimesShared += amount;
    shape->numTimesShared += amount;
  };



public:

  // a gamma distribution always has dimensionalit of 1 (but we do
  // a vector of Gammas, sort of like a Diagonal component Gaussian.
  GammaComponent(const int dim) : Component(dim) { }

  virtual ~GammaComponent() { }

  static double varianceFloor() { return _varianceFloor; }
  static double setVarianceFloor(const double floor);

  static bool cloneShareScale;
  static bool cloneShareShape;

  //////////////////////////////////////////////
  // read/write basic parameters
  void read(iDataStreamFile& is);
  void write(oDataStreamFile& os);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////
  // create a copy of self, with entirel or partially new parameters
  // (so clone might or might not share, depends on the object), and
  // with slightly (and randomly) perturbed values.
  Component* noisyClone();

  /////////////////////////////////////////////////
  // create a copy of self, with entirely new parameters with
  // identical values; NOTHING is shared
  Component* identicalIndependentClone();

  //////////////////////////////////
  // set all current parameters to valid but random values
  void makeRandom();
  // set all current parameters to valid but "uniform" values 
  // (for Gammas this means shape = 1, and scale = 2.
  void makeUniform();
  unsigned totalNumberParameters() { return 2; }
  //////////////////////////////////

  void recursivelyClearUsedBit() { 
    emClearUsedBit(); 
    scale->recursivelyClearUsedBit();
    shape->recursivelyClearUsedBit();
  } 
  void recursivelySetUsedBit() { 
    emSetUsedBit(); 
    scale->recursivelySetUsedBit();
    shape->recursivelySetUsedBit();
  }

  //////////////////////////////////
  // probability evaluation
  logpr log_p(const float *const x,     // real-valued scoring obs at time t
		      const Data32* const base, // ptr to base obs at time t
		      const int stride);    // stride

  // TODO: finish this next function, returning the score of the mode of the distribution.
  // virtual logpr maxValue() { error("Gamma max value not implemented\n"); }
  //////////////////////////////////


  //////////////////////////////////
  // Full Baum-Welch EM training  //
  //////////////////////////////////
  void emStartIteration();
  void emIncrement(logpr prob,
		   const float*f,
		   const Data32* const base,
		   const int stride);
  void emEndIteration();
  void emSwapCurAndNew();

  // parallel training
  void emStoreObjectsAccumulators(oDataStreamFile& ofile,
				  bool writeLogVals = true,
				  bool writeZeros = false);
  void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile);
  void emZeroOutObjectsAccumulators();
  void emLoadObjectsAccumulators(iDataStreamFile& ifile);
  void emAccumulateObjectsAccumulators(iDataStreamFile& ifile);
  const string typeName() { return "Gamma Distribution"; }
  //////////////////////////////////

  //////////////////////////////////
  // Sample Generation            //
  //////////////////////////////////
  virtual void sampleGenerate(float *sample,
			      const Data32* const base,
			      const int stride);
  //////////////////////////////////

};


#endif
